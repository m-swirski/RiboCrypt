# RiboCrypt 1.15.5

## Experimental feature unlock (`experimental_feature_password`)

* New argument `RiboCrypt_app(experimental_unlock = NULL)`: a shared unlock code for
  private / experimental features. No user accounts: collaborators type the code in
  the masked "Experimental features (code)" box in the browser or observatory
  settings, and the features unlock for their own browser session only (a reload
  locks again; nothing is stored in URLs, bookmarks or cookies).
* The code can also come from `options(RiboCrypt.experimental_unlock_file = path)`
  (recommended, together with `experimental_unlock_file(path)`; file mode 600),
  `options(RiboCrypt.experimental_unlock = value)` or the environment variables
  `RIBOCRYPT_EXPERIMENTAL_UNLOCK_FILE` / `RIBOCRYPT_EXPERIMENTAL_UNLOCK`.
  Only a salted HMAC-SHA256 of the code is kept in memory; plain option and
  environment values are cleared after reading. The code never reaches the client
  (HTML, outputs, share URLs, cache keys or error messages).
* Wrong codes are rate limited (2 s cool-down, at most 5 wrong codes per minute per
  session, 50 per minute per process). Feedback is only "Unlocked" / "Wrong code".
* Generic gate for later features: `experimental_feature_enabled(session, "<feature>")`.

## Behaviour change: explainability mode needs the unlock

* The explainability ("expectancy") mode is now an experimental feature. The
  explain store option (`RiboCrypt.explain_store` / `RIBOCRYPT_EXPLAIN_STORE`) alone
  is no longer enough: the explain settings (status, expected E, pLDDT, UniProt,
  translon and Paralloc lanes) are rendered and computed only for sessions unlocked
  with the code. Without `experimental_unlock` configured (default) the mode stays
  hidden, also for local use; for local testing pass any value and type it in.
* `explain_mode=TRUE` in a shared URL is ignored for locked sessions and is only
  written to share URLs from unlocked sessions.
