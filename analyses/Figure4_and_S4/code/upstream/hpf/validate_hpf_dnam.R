#!/usr/bin/env Rscript
# Historical RDS files are read only here, after independent regeneration.
args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=3L) stop("Usage: Rscript validate_hpf_dnam.R REBUILT_DIR HISTORICAL_DIR REPORT_DIR")
rebuilt_dir <- normalizePath(args[1],mustWork=TRUE)
historical_dir <- normalizePath(args[2],mustWork=TRUE)
report_dir <- args[3]
dir.create(report_dir,recursive=TRUE,showWarnings=FALSE)
files <- c("HPF_probe_retention_table_full.rds", "HPF_DNAm_retention_decision_trials.rds",
           "HPF_DNAm_FINAL_CORE_TRIALS.rds")
checks <- list(); attributes_log <- list(); metadata_log <- list(); count <- 0L
allowed_metadata <- c("HPF_DNAm_retention_decision_trials.rds/source/note",
                     "HPF_DNAm_FINAL_CORE_TRIALS.rds/source/probe_table")
record <- function(path,n,exact,error,passed,kind) {
  count <<- count+1L
  checks[[count]] <<- data.frame(path=path,n=n,identical_values=exact,
    maximum_absolute_error=error,passed=passed,kind=kind,stringsAsFactors=FALSE)
}
compare <- function(actual,expected,path) {
  if(path %in% allowed_metadata) {
    metadata_log[[length(metadata_log)+1L]] <<- data.frame(path=path,
      historical=paste(expected,collapse=" | "),rebuilt=paste(actual,collapse=" | "))
    return(invisible(NULL))
  }
  at <- attributes(actual); et <- attributes(expected)
  if(!identical(at,et)) attributes_log[[length(attributes_log)+1L]] <<-
    data.frame(path=path,detail=paste(all.equal(at,et),collapse="; "))
  if(is.list(expected)) {
    names_ok <- identical(names(actual),names(expected)) && length(actual)==length(expected)
    record(paste0(path,"/<structure>"),length(expected),names_ok,NA_real_,names_ok,"structure")
    if(!names_ok) return(invisible(NULL))
    for(nm in names(expected)) compare(actual[[nm]],expected[[nm]],paste(path,nm,sep="/"))
  } else if(is.numeric(expected) && is.numeric(actual)) {
    compatible <- length(actual)==length(expected) && identical(is.na(actual),is.na(expected)) &&
      identical(is.finite(actual),is.finite(expected))
    err <- if(compatible && any(is.finite(expected))) max(abs(actual[is.finite(expected)]-expected[is.finite(expected)])) else if(compatible) 0 else Inf
    exact <- identical(unname(actual),unname(expected))
    record(path,length(expected),exact,err,compatible && err<=1e-10,"numeric")
  } else {
    exact <- identical(unname(actual),unname(expected))
    record(path,length(expected),exact,NA_real_,exact,"categorical")
  }
}
for(file in files) {
  message("Comparing ",file)
  compare(readRDS(file.path(rebuilt_dir,file)),readRDS(file.path(historical_dir,file)),file)
}
out <- do.call(rbind,checks)
write.table(out,file.path(report_dir,"all_fields_comparison.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
if(length(attributes_log)) write.table(do.call(rbind,attributes_log),file.path(report_dir,"attribute_differences.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
if(length(metadata_log)) write.table(do.call(rbind,metadata_log),file.path(report_dir,"documented_metadata_changes.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
summary <- do.call(rbind,lapply(files,function(file) {
  a <- out[startsWith(out$path,paste0(file,"/")),]
  data.frame(file=file,n_checks=nrow(a),passed=sum(a$passed),failed=sum(!a$passed),
    maximum_absolute_error=max(a$maximum_absolute_error,na.rm=TRUE))
}))
write.table(summary,file.path(report_dir,"comparison_summary.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
print(summary,row.names=FALSE)
if(any(!out$passed)) {print(out[!out$passed,]);stop("Analytical comparison failed; see the report.")}
message("All analytical values, IDs and ordering match within absolute tolerance 1e-10.")
