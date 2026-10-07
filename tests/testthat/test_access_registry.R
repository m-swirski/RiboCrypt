access_test_registry <- function(path = ":memory:") {
  skip_if_not_installed("DBI")
  skip_if_not_installed("RSQLite")
  con <- ribocrypt_access_db(path)
  ribocrypt_access_workspace(con, "a")
  ribocrypt_access_workspace(con, "b")
  ribocrypt_access_member(con, "issuer", "alice", "a")
  ribocrypt_access_member(con, "issuer", "bob", "b")
  for (id in c("public", "private-a", "private-b")) {
    directory <- file.path(tempdir(), paste0("ribocrypt-access-", id))
    dir.create(directory, showWarnings = FALSE)
    ribocrypt_access_dataset(con, id, "same-name", directory, public = id == "public")
  }
  ribocrypt_access_grant(con, "a", "private-a")
  ribocrypt_access_grant(con, "b", "private-b", download = TRUE)
  con
}

test_that("registry separates users, issuers and capabilities", {
  con <- access_test_registry()
  on.exit(DBI::dbDisconnect(con))
  alice <- list(issuer = "issuer", subject = "alice")
  bob <- list(issuer = "issuer", subject = "bob")
  expect_identical(ribocrypt_access_catalog(con)$id, "public")
  expect_identical(ribocrypt_access_catalog(con, alice)$id, c("private-a", "public"))
  expect_identical(ribocrypt_access_catalog(con, bob)$id, c("private-b", "public"))
  expect_identical(ribocrypt_access_catalog(con, alice, "download")$id, "public")
  expect_identical(ribocrypt_access_catalog(con, bob, "download")$id, c("private-b", "public"))
  expect_identical(ribocrypt_access_catalog(con, list(issuer = "other", subject = "alice"))$id, "public")
  expect_identical(ribocrypt_access_catalog(con, list(issuer = "issuer", subject = "new-account"))$id, "public")
  expect_error(ribocrypt_access_resolve(con, "private-b", alice), "unavailable or access denied")
  expect_error(ribocrypt_access_resolve(con, "missing", alice), "unavailable or access denied")
  expect_error(ribocrypt_access_resolve(con, "same-name", bob), "unavailable or access denied")
  expect_error(ribocrypt_access_resolve(con, "private-a", alice, "download"), "access denied")
  expect_identical(ribocrypt_access_resolve(con, "private-b", bob)$experiment, "same-name")
  expect_false(identical(ribocrypt_access_resolve(con, "private-a", alice)$directory,
                         ribocrypt_access_resolve(con, "private-b", bob)$directory))
})

test_that("revocation is checked on every resolution and membership is required", {
  con <- access_test_registry()
  on.exit(DBI::dbDisconnect(con))
  alice <- list(issuer = "issuer", subject = "alice")
  expect_equal(nrow(ribocrypt_access_resolve(con, "private-a", alice)), 1)
  ribocrypt_access_grant(con, "a", "private-a", enabled = FALSE)
  expect_error(ribocrypt_access_resolve(con, "private-a", alice), "access denied")
  ribocrypt_access_grant(con, "a", "private-a", download = TRUE)
  expect_equal(nrow(ribocrypt_access_resolve(con, "private-a", alice, "download")), 1)
  ribocrypt_access_member(con, "issuer", "alice", "a", enabled = FALSE)
  expect_error(ribocrypt_access_resolve(con, "private-a", alice), "access denied")
})

test_that("persistent registry survives reopening without changing mappings", {
  path <- tempfile(fileext = ".sqlite")
  on.exit(unlink(path))
  con <- access_test_registry(path)
  DBI::dbDisconnect(con)
  con <- ribocrypt_access_db(path)
  on.exit(DBI::dbDisconnect(con), add = TRUE)
  alice <- list(issuer = "issuer", subject = "alice")
  expect_identical(ribocrypt_access_catalog(con, alice)$id, c("private-a", "public"))
  expect_error(ribocrypt_access_dataset(con, "private-a", "different", tempdir()), "UNIQUE")
  expect_identical(ribocrypt_access_resolve(con, "private-a", alice)$experiment, "same-name")
})

test_that("registry rejects malformed inputs and unknown references", {
  skip_if_not_installed("DBI")
  skip_if_not_installed("RSQLite")
  con <- ribocrypt_access_db(":memory:")
  on.exit(DBI::dbDisconnect(con))
  expect_error(ribocrypt_access_workspace(con, ""), "nonempty")
  expect_error(ribocrypt_access_member(con, "issuer", "user", "missing"), "FOREIGN KEY")
  expect_error(ribocrypt_access_dataset(con, "x", "../private", tempdir()), "not a path")
  expect_error(ribocrypt_access_dataset(con, "x", "..\\private", tempdir()), "not a path")
  expect_error(ribocrypt_access_dataset(con, "x", "valid", tempfile()), "does not exist")
  expect_error(ribocrypt_access_dataset(con, "x", "valid", tempdir(), NA), "TRUE or FALSE")
  expect_error(ribocrypt_access_workspace(con, c("a", "b")), "nonempty")
  ribocrypt_access_workspace(con, "quote' OR 1=1 --")
  expect_equal(nrow(DBI::dbReadTable(con, "workspaces")), 1)
  expect_error(ribocrypt_access_catalog(con, list(issuer = "issuer")), "subject")
})

test_that("another connection's revocation is immediately visible", {
  path <- tempfile(fileext = ".sqlite")
  on.exit(unlink(path))
  administrator <- access_test_registry(path)
  session <- ribocrypt_access_db(path)
  on.exit(DBI::dbDisconnect(administrator), add = TRUE)
  on.exit(DBI::dbDisconnect(session), add = TRUE)
  alice <- list(issuer = "issuer", subject = "alice")
  expect_equal(nrow(ribocrypt_access_resolve(session, "private-a", alice)), 1)
  ribocrypt_access_grant(administrator, "a", "private-a", enabled = FALSE)
  expect_error(ribocrypt_access_resolve(session, "private-a", alice), "access denied")
  expect_error(ribocrypt_access_resolve(session, "private-a' OR 1=1 --", alice), "access denied")
})

test_that("gateway verification rejects forged identity and fixes the issuer", {
  secret <- paste(rep("a", 64), collapse = "")
  request <- list(HTTP_RIBOCRYPT_GATEWAY_SECRET = secret, HTTP_RIBOCRYPT_SUBJECT = "user-123")
  expect_identical(ribocrypt_access_identity(request, "fixed-issuer", secret),
                   list(issuer = "fixed-issuer", subject = "user-123"))
  expect_error(ribocrypt_access_identity(list(HTTP_RIBOCRYPT_SUBJECT = "admin"), "issuer", secret), "Untrusted")
  request$HTTP_RIBOCRYPT_GATEWAY_SECRET <- paste0("b", substring(secret, 2))
  expect_error(ribocrypt_access_identity(request, "issuer", secret), "Untrusted")
  request$HTTP_RIBOCRYPT_GATEWAY_SECRET <- secret
  request$HTTP_RIBOCRYPT_SUBJECT <- ""
  expect_null(ribocrypt_access_identity(request, "issuer", secret))
  request$HTTP_RIBOCRYPT_SUBJECT <- "user\nadmin"
  expect_error(ribocrypt_access_identity(request, "issuer", secret), "control characters")
  expect_error(ribocrypt_access_identity(request, "issuer", "short"), "32 bytes")
})
