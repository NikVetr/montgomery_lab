x <- MotrpacBicQC::tissue_abbr

x <- x[!grepl(names(x), pattern = "-")]
xn <- unique(x)
x <- x[match(xn, x)]

nogo <- c("VENACV", "PLASMA", "OVARY", "TESTES")

x <- x[!(x %in% nogo)]


cat(paste0("\\textcolor[HTML]{", gsub("#", "", MotrpacBicQC::tissue_cols[x]), "}{", x, "} & ", stringr::str_to_title(names(x)), "\\\\\n"), sep = "")


