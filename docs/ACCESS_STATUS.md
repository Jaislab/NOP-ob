# Original cloud environment: access status

Last checked November 2, 2026, approximately 09:55 CET (08:55 UTC).

| Check | Observed outcome |
|---|---|
| Current desired/observed configuration revision | Both revision 4 |
| Configured internet mode | Unrestricted |
| Reported network-enforcement state | Unknown |
| www.ncbi.nlm.nih.gov via inherited HTTPS proxy | CONNECT rejected, HTTP 403 Forbidden |
| ftp.ncbi.nlm.nih.gov via inherited HTTPS proxy | CONNECT rejected, HTTP 403 Forbidden |
| github.com via the same proxy | HTTP 200 |
| registry.npmjs.org via the same proxy | HTTP 200 |
| Retained diagnostic records | None available |

Requests included the GSE225224 GEO text record, its supplementary directory
listing and its family SOFT archive. The failures occurred at the proxy connection
stage, so they do not establish whether any listed dataset exists or is public.
The cause of the discrepancy between saved configuration and actual proxy
behavior was not established. No credentials, raw proxy logs or private runtime
identifiers are included in this repository.

The next executor should check actual official GEO access from its own authorized
environment and proceed with acquisition when it works. The original blocker
does not need to be reproduced to perform the biological study.
