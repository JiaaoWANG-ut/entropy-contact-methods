#!/usr/bin/env python3
"""Inspect the Gmail Sent folder over IMAP, tolerating the flaky egress path.

Connections are made per resolved IP with the correct TLS SNI hostname, because
only part of Google's address pool is routable from the cluster.
"""

from __future__ import annotations

import imaplib
import os
import socket
import ssl
import sys
import time

HOST = "imap.gmail.com"
PORT = 993
USER = "123jiaao@gmail.com"


class PinnedIMAP(imaplib.IMAP4_SSL):
    def __init__(self, ip: str, timeout: int = 20):
        self._ip = ip
        super().__init__(host=HOST, port=PORT, timeout=timeout,
                         ssl_context=ssl.create_default_context())

    def _create_socket(self, timeout):
        sock = socket.create_connection((self._ip, PORT), timeout=timeout)
        return self.ssl_context.wrap_socket(sock, server_hostname=HOST)

    def open(self, host="", port=PORT, timeout=None):
        self.host, self.port = HOST, PORT
        self.sock = self._create_socket(timeout or 20)
        self.file = self.sock.makefile("rb")


def ips():
    return sorted({i[4][0] for i in socket.getaddrinfo(
        HOST, PORT, socket.AF_INET, socket.SOCK_STREAM)})


def connect(password: str, deadline_s: int = 600):
    end = time.time() + deadline_s
    while time.time() < end:
        for ip in ips():
            try:
                imap = PinnedIMAP(ip)
                imap.login(USER, password)
                print(f"IMAP connected via {ip}")
                return imap
            except Exception as exc:
                print(f"  {ip}: {type(exc).__name__}", flush=True)
        time.sleep(5)
    raise SystemExit("IMAP unreachable")


def main():
    password = os.environ["GMAIL_APP_PASSWORD"]
    query = sys.argv[1] if len(sys.argv) > 1 else "entropy analysis"
    imap = connect(password)
    try:
        imap.select('"[Gmail]/Sent Mail"', readonly=True)
        typ, data = imap.search(None, "SUBJECT", f'"{query}"')
        nums = data[0].split()
        print(f"Matching sent messages: {len(nums)}")
        for num in nums[-5:]:
            typ, payload = imap.fetch(
                num, "(BODY[HEADER.FIELDS (DATE TO SUBJECT)])")
            print(payload[0][1].decode(errors="replace").strip())
            print("-" * 60)
    finally:
        try:
            imap.logout()
        except Exception:
            pass


if __name__ == "__main__":
    main()
