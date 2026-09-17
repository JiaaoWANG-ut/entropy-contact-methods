#!/usr/bin/env python3
"""
Mail an entropy analysis directory (figures + tables) to the group.

Only part of the Gmail address pool is routable from this cluster and the link
can drop mid-DATA, so each resolved IP is tried in turn and, before any retry,
the Sent folder is checked for the message's Message-Id to avoid delivering a
duplicate.
"""

from __future__ import annotations

import argparse
import os
import smtplib
import socket
import ssl
import time
import zipfile
from email.message import EmailMessage
from email.utils import formatdate, make_msgid

SMTP_HOST = "smtp.gmail.com"
SMTP_PORTS = (465, 587)
SENDER = "123jiaao@gmail.com"
APP_PASSWORD = os.environ.get("GMAIL_APP_PASSWORD", "")
RECIPIENTS = [
    "lq.wang@utoronto.ca",
    "xue.yao@utoronto.ca",
    "yimin.wu@uwaterloo.ca",
    "wangjiaao0720@utexas.edu",
]

ROOT = os.path.dirname(os.path.abspath(__file__))


def resolved_ips(host: str, port: int):
    infos = socket.getaddrinfo(host, port, socket.AF_INET, socket.SOCK_STREAM)
    seen, ips = set(), []
    for info in infos:
        ip = info[4][0]
        if ip not in seen:
            seen.add(ip)
            ips.append(ip)
    return ips


# Addresses that have carried a successful session from this cluster; DNS does
# not always hand them out, so keep them in the pool.
KNOWN_GOOD_IPS = ("64.233.188.108", "64.233.188.109", "64.233.187.108",
                  "64.233.187.109")


def candidate_ips(port: int):
    """Accumulate IPs from both the public name and its CNAME target."""
    ips = list(KNOWN_GOOD_IPS)
    for host in (SMTP_HOST, "gmail-smtp-msa.l.google.com"):
        try:
            for ip in resolved_ips(host, port):
                if ip not in ips:
                    ips.append(ip)
        except socket.gaierror:
            pass
    return ips


def connect(deadline_s: int = 2700, timeout: int = 60):
    """Return a logged-in SMTP connection, retrying until the deadline."""
    context = ssl.create_default_context()
    errors = []
    end = time.time() + deadline_s
    round_no = 0
    while time.time() < end:
        round_no += 1
        for port in SMTP_PORTS:
            for ip in candidate_ips(port):
                try:
                    sock = socket.create_connection((ip, port), timeout=timeout)
                    if port == 465:
                        sock = context.wrap_socket(sock, server_hostname=SMTP_HOST)
                    server = (smtplib.SMTP_SSL if port == 465 else smtplib.SMTP)(
                        host=SMTP_HOST, timeout=timeout)
                    server.sock = sock
                    code, _ = server.getreply()
                    if code != 220:
                        raise smtplib.SMTPConnectError(code, "no greeting")
                    server.ehlo()
                    if port != 465:
                        server.starttls(context=context)
                        server.ehlo()
                    server.login(SENDER, APP_PASSWORD)
                    print(f"SMTP connected via {ip}:{port}")
                    return server
                except smtplib.SMTPAuthenticationError:
                    raise
                except Exception as exc:  # unroutable IP, TLS reset, timeout
                    errors.append(f"{ip}:{port} {type(exc).__name__}: {exc}")
        print(f"  round {round_no}: no route yet, retrying", flush=True)
        time.sleep(5)
    raise SystemExit("SMTP unreachable:\n  " + "\n  ".join(errors[-12:]))


def already_sent(message_id: str) -> bool:
    """True if message_id is present in the Gmail Sent folder."""
    try:
        import check_sent

        imap = check_sent.connect(APP_PASSWORD, deadline_s=240)
    except SystemExit:
        print("  could not verify via IMAP")
        return False
    try:
        imap.select('"[Gmail]/Sent Mail"', readonly=True)
        typ, data = imap.search(None, "HEADER", "Message-ID", message_id)
        return bool(data and data[0].split())
    except Exception as exc:
        print(f"  IMAP check failed: {type(exc).__name__}: {exc}")
        return False
    finally:
        try:
            imap.logout()
        except Exception:
            pass


def send_with_retry(msg, attempts: int = 8):
    message_id = msg["Message-ID"]
    for attempt in range(1, attempts + 1):
        server = connect()
        try:
            server.send_message(msg)
            print(f"Message accepted by Gmail (attempt {attempt})")
            return
        except Exception as exc:
            print(f"  attempt {attempt} failed: {type(exc).__name__}: {exc}",
                  flush=True)
        finally:
            try:
                server.quit()
            except Exception:
                pass
        print("  checking whether it landed anyway ...", flush=True)
        if already_sent(message_id):
            print("Message is in the Sent folder; not resending.")
            return
        time.sleep(5)
    raise SystemExit(f"Could not deliver the message after {attempts} attempts.")


def collect_files(out_dir: str):
    figures = sorted(os.path.join(out_dir, f) for f in os.listdir(out_dir)
                     if f.endswith((".png", ".svg")))
    tables = sorted(os.path.join(out_dir, f) for f in os.listdir(out_dir)
                    if f.endswith((".csv", ".xlsx")))
    return figures, tables


def make_zip(out_dir: str, zip_name: str, files):
    path = os.path.join(out_dir, zip_name)
    stem = os.path.splitext(zip_name)[0]
    with zipfile.ZipFile(path, "w", zipfile.ZIP_DEFLATED) as zf:
        for src in files:
            zf.write(src, os.path.join(stem, os.path.basename(src)))
    return path


def build_message(subject: str, body: str, attachments):
    msg = EmailMessage()
    msg["From"] = f"Jiaao Wang <{SENDER}>"
    msg["To"] = ", ".join(RECIPIENTS)
    msg["Subject"] = subject
    msg["Date"] = formatdate(localtime=True)
    msg["Message-ID"] = make_msgid(domain="gmail.com")
    msg.set_content(body)
    for path in attachments:
        with open(path, "rb") as fh:
            data = fh.read()
        name = os.path.basename(path)
        if name.endswith(".png"):
            maintype, subtype = "image", "png"
        elif name.endswith(".zip"):
            maintype, subtype = "application", "zip"
        elif name.endswith(".xlsx"):
            maintype, subtype = ("application",
                                 "vnd.openxmlformats-officedocument."
                                 "spreadsheetml.sheet")
        else:
            maintype, subtype = "text", "plain"
        msg.add_attachment(data, maintype=maintype, subtype=subtype,
                           filename=name)
    return msg


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out-dir", default=os.path.join(ROOT, "entropy_analysis"))
    parser.add_argument("--zip-name", default="FIG4_entropy_analysis.zip")
    parser.add_argument("--subject", default=(
        "FIG4 NiO-6/NiO-8 entropy analysis: S_atom and S_config vs time "
        "(1500-3000 K)"))
    parser.add_argument("--body-file", default=None)
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()

    if not APP_PASSWORD:
        raise SystemExit("Set GMAIL_APP_PASSWORD in the environment.")

    out_dir = os.path.abspath(args.out_dir)
    body_file = args.body_file or os.path.join(out_dir, "email_body.txt")

    figures, tables = collect_files(out_dir)
    archive = make_zip(out_dir, args.zip_name, figures + tables)
    attachments = ([f for f in figures if f.endswith(".png")]
                   + [t for t in tables if t.endswith(".xlsx")]
                   + [archive])

    with open(body_file, encoding="utf-8") as fh:
        body = fh.read()

    msg = build_message(args.subject, body, attachments)
    total_mb = sum(os.path.getsize(p) for p in attachments) / 1e6
    print(f"Source: {out_dir}")
    print(f"Attachments: {len(attachments)} files, {total_mb:.1f} MB")
    for path in attachments:
        print(f"  {os.path.basename(path)}  {os.path.getsize(path)/1e6:.2f} MB")

    if args.dry_run:
        print("Dry run: not sending.")
        return

    send_with_retry(msg)
    print("Sent to: " + ", ".join(RECIPIENTS))


if __name__ == "__main__":
    main()
