# Local example: Rscript run.R S4 baseline CC 1
# Production: one setting, condition and method per task; never change the solver.
source("generate_data.R")
source("methods.R")

# Compute diagnostics after fitting. No oracle quantity is passed to fit_itr().
evaluate_fit <- function(sample, method, fit) {
  d <- paper_setting(sample$setting)
  test <- sample$test
  treated <- fit$theta[1] + fit$theta[2]*test$x1 > 0
  tau <- d$tau(test$x1,test$x2)
  errors <- list(q_RMSE=NA_real_, q_best_constant_RMSE=NA_real_,
                 q_oracle=NA_real_, q_fitted=NA_real_)
  if (method == "Fusion") {
    X <- sample$target
    sigma1 <- conditions$sigma1[match(sample$condition,conditions$condition)]
    density <- function(sd)
      dnorm(X$x1,sd=sd)*dnorm(X$x2,sd=sd)/(pnorm(3/sd)-pnorm(-3/sd))^2
    w <- density(1)/density(sigma1)
    trial_noise <- w*d$eta1^2/nrow(sample$trial)
    target_noise <- (d$eta0(X$x1)^2 + .75*d$c(X$x1,X$x2)^2)/nrow(X)
    q_star <- trial_noise/(trial_noise+target_noise)
    stopifnot(length(fit$q_hat)==nrow(X), all(is.finite(fit$q_hat)))
    errors <- list(q_RMSE=sqrt(mean((fit$q_hat-q_star)^2)),
                   q_best_constant_RMSE=sqrt(mean((q_star-mean(q_star))^2)),
                   q_oracle=mean(q_star), q_fitted=mean(fit$q_hat))
  }
  data.frame(method,variant=sample$condition,setting=sample$setting,
             rep=sample$rep,seed=sample$seed,solver="DC",
             dc_starts=fit$dc_starts,bfgs_calls=fit$bfgs_calls,
             regret=population_regret(fit$theta,sample$setting),
             test_regret=mean(((tau>0)-treated)*tau),
             intercept=fit$theta[1],slope=fit$theta[2],
             rank_deficient_fits=fit$rank_deficient_fits, errors)
}

out_dir <- Sys.getenv("ITR_RESULTS","results")
result_path <- function(setting,condition,method,n)
  file.path(out_dir,sprintf("%s_%s_%s_r%d.csv",setting,condition,method,n))

run_cell <- function(setting, condition, method, n) {
  stopifnot(setting %in% c("S1","S2","S3","S4"), condition %in% conditions$condition,
            method %in% method_names, n %in% 1:25)
  path <- result_path(setting,condition,method,n)
  if (file.exists(path)) stop("Output already exists: ",path,". Choose another ITR_RESULTS folder.")
  dir.create(out_dir,recursive=TRUE,showWarnings=FALSE)
  results <- lapply(seq_len(n),function(rep) {
    sample <- generate_data(setting,condition,rep)
    fit <- fit_itr(sample,method)
    evaluate_fit(sample,method,fit)
  })
  rows <- do.call(rbind,results)
  stopifnot(all(is.finite(rows$regret)), all(rows$regret>=0))
  tmp <- tempfile(tmpdir=out_dir)
  write.csv(rows,tmp,row.names=FALSE); stopifnot(file.rename(tmp,path))
  message(path)
  invisible(rows)
}

# A single array contains only S4: six conditions, nine methods, 25 paired seeds.
tasks <- expand.grid(method=method_names,condition=conditions$condition,stringsAsFactors=FALSE)
collect_s4 <- function() {
  paths <- mapply(result_path,"S4",tasks$condition,tasks$method,25L)
  if (any(!file.exists(paths))) stop("Missing cells: ",paste(which(!file.exists(paths)),collapse=", "))
  rows <- do.call(rbind,lapply(paths,read.csv))
  stopifnot(nrow(rows)==1350L, all(rows$setting=="S4"), all(rows$seed==9730000L+rows$rep),
            setequal(rows$rep,1:25),setequal(rows$method,method_names),
            setequal(rows$variant,conditions$condition),all(rows$solver=="DC"),
            all(is.finite(rows$regret)),all(rows$regret>=0),all(rows$dc_starts>=16),
            all(rows$bfgs_calls>=rows$dc_starts),all(table(rows$method,rows$variant,rows$rep)==1))
  put <- function(x,name) write.csv(x,file.path(out_dir,paste0("s4_r25_",name,".csv")),row.names=FALSE)
  put(rows,"replications")
  all_methods <- do.call(rbind,lapply(split(rows,list(rows$method,rows$variant)),function(r)
    data.frame(method=r$method[1],variant=r$variant[1],n=nrow(r),mean_regret=mean(r$regret),
               sd_regret=sd(r$regret),median_regret=median(r$regret),max_regret=max(r$regret))))
  put(all_methods,"all_methods_summary")
  paired <- do.call(rbind,lapply(conditions$condition,function(condition) {
    r <- rows[rows$variant==condition,]
    regret <- function(method) {a <- r[r$method==method,]; a$regret[order(a$rep)]}
    gap <- regret("CC")-regret("IP")
    data.frame(variant=condition,WY=mean(regret("WY")),IP=mean(regret("IP")),
               CC=mean(regret("CC")),Fusion=mean(regret("Fusion")),
               CC_minus_IP=mean(gap),MCSE=sd(gap)/sqrt(25),median=median(gap),
               CC_wins=sum(gap<0),ties=sum(gap==0),
               largest_gap=max(gap),max_leave_one_out=max((sum(gap)-gap)/24))
  }))
  put(paired,"paired_summary"); print(paired,row.names=FALSE)
}

compare_s2_s4 <- function(previous_csv) {
  old <- read.csv(previous_csv)
  new <- read.csv(file.path(out_dir,"s4_r25_replications.csv"))
  gap <- function(x,setting) {
    x <- x[x$setting==setting & x$method %in% c("IP","CC"),]
    w <- reshape(x[c("variant","rep","seed","method","regret")],
                 idvar=c("variant","rep","seed"),timevar="method",direction="wide")
    stopifnot(nrow(w)==150L,!anyNA(w))
    transform(w,delta=regret.CC-regret.IP)[c("variant","rep","seed","delta")]
  }
  paired <- merge(gap(old,"S2a_quadratic_c"),gap(new,"S4"),
                  by=c("variant","rep","seed"),suffixes=c("_S2","_S4"))
  stopifnot(nrow(paired)==150L,all(table(paired$variant)==25L))
  print(do.call(rbind,lapply(split(paired,paired$variant),function(r){
    change <- r$delta_S4-r$delta_S2
    data.frame(condition=r$variant[1],S2=mean(r$delta_S2),S4=mean(r$delta_S4),
               change=mean(change),paired_MCSE=sd(change)/sqrt(nrow(r)))
  })),row.names=FALSE)
}

verify_paper_fit <- function(csv, setting, method, rep) {
  paper_ids <- c(S1="S1_IP",S2="S2a_quadratic_c",S3="S3_Fusion")
  stopifnot(setting %in% names(paper_ids), method %in% method_names, rep %in% 1:25)
  ref <- read.csv(csv)
  ref <- ref[ref$setting==paper_ids[[setting]] & ref$variant=="baseline" &
               ref$method==method & ref$rep==rep,]
  sample <- generate_data(setting,"baseline",rep)
  stopifnot(nrow(ref)==1L, ref$seed==sample$seed)
  new <- evaluate_fit(sample,method,fit_itr(sample,method))
  fields <- c("regret","test_regret","intercept","slope","dc_starts","bfgs_calls",
              "rank_deficient_fits",if(method=="Fusion") "q_RMSE")
  check <- data.frame(metric=fields, archived=as.numeric(unlist(ref[fields])),
                      rerun=as.numeric(unlist(new[fields])))
  check$difference <- check$rerun-check$archived
  print(check,row.names=FALSE,digits=12)
  if (any(abs(check$difference)>1e-9))
    warning("Numerical differences from the archive. Check R, packages, BLAS and solver diagnostics.")
}

if (sys.nframe()==0L) {
  args <- commandArgs(TRUE)
  if (!length(args)) stop("Example: Rscript run.R S4 baseline CC 1")
  if (args[1]=="manifest") {
    print(data.frame(task=seq_len(nrow(tasks)),tasks),row.names=FALSE)
  } else if (args[1]=="task") {
    id <- as.integer(args[2]); stopifnot(id %in% seq_len(nrow(tasks)))
    run_cell("S4",tasks$condition[id],tasks$method[id],25L)
  } else if (args[1]=="collect") {
    collect_s4()
  } else if (args[1]=="compare") {
    compare_s2_s4(args[2])
  } else if (args[1]=="verify") {
    stopifnot(length(args)==5L)
    verify_paper_fit(args[2],args[3],args[4],as.integer(args[5]))
  } else {
    stopifnot(length(args)==4L)
    requested <- if(args[3]=="all") method_names else args[3]
    for(method in requested) run_cell(args[1],args[2],method,as.integer(args[4]))
  }
}
