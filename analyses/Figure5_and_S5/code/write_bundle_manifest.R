# Inventory the delivered files, excluding this manifest's own checksum.
files <- sort(list.files(".",recursive=TRUE,full.names=FALSE,all.files=FALSE))
files <- files[!file.info(files)$isdir & files!="provenance/bundle_checksums.tsv"]
write.table(data.frame(File=files,Bytes=file.info(files)$size,
  MD5=unname(tools::md5sum(files))),
  "provenance/bundle_checksums.tsv",sep="\t",quote=FALSE,row.names=FALSE)
