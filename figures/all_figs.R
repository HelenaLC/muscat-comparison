scripts <- list.files("figures", "*.R$", full.names = TRUE)
scripts <- scripts[!grepl("all_figs.R", scripts)]
for (s in scripts) tryCatch(source(s), error = function(e) message(e))
