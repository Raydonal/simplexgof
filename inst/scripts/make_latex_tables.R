outdir<-"/home/claude/paper_outputs"

# ── Table: Un characteristic measures ───────────────────────────
meas<-read.csv(file.path(outdir,"table_Un_measures.csv"))
mu_tex<-c(low="$(0.019,0.147)$",mid="$(0.205,0.886)$",high="$(0.903,0.995)$")
tex<-c("\\begin{table}[htb!]","\\centering",
"\\caption{Estimated characteristic measures of the $U_n$ distribution, $n=40$, based on 1000 Monte Carlo replications.}\\label{tab:Un_measures}",
"\\begin{tabular}{cc|cccc}\\hline",
"$\\sigma^2$ & $\\mu$ range & Mean & Variance & Kurtosis & Skewness \\\\\\hline")
for(i in 1:nrow(meas)){r<-meas[i,]
  tex<-c(tex,sprintf("%.1f & %s & %.3f & %.3f & %.3f & %.3f \\\\",
    r$sigma2,mu_tex[r$mu_range],r$Mean,r$Variance,r$Kurtosis,r$Skewness))
  if(i==3)tex<-c(tex,"\\hline")}
tex<-c(tex,"\\hline","\\end{tabular}","\\end{table}")
writeLines(tex,file.path(outdir,"table_Un_measures.tex"))

# ── Application param tables (LaTeX) ────────────────────────────
for(app in c("ammonia","pbsc")){
  pp<-read.csv(file.path(outdir,paste0("table_",app,"_params.csv")))
  gg<-read.csv(file.path(outdir,paste0("table_",app,"_gof.csv")))
  cap<-if(app=="ammonia")"ammonia oxidation" else "PBSC transplant"
  texp<-c("\\begin{table}[htb!]","\\centering",
    sprintf("\\caption{Maximum likelihood estimates for the %s data.}\\label{tab:%s_par}",cap,app),
    "\\begin{tabular}{l l rrr}\\hline",
    "Parameter & Sub-model & Estimate & Std.\\ Error & $p$-value \\\\\\hline")
  for(i in 1:nrow(pp)){r<-pp[i,]
    pname<-gsub("beta","\\\\beta_",gsub("gamma","\\\\gamma_",r$Parameter))
    texp<-c(texp,sprintf("$%s$ & %s & %.4f & %.4f & %s \\\\",
      pname,r$Sub_model,r$Estimate,r$Std_Error,r$p_value))}
  texp<-c(texp,"\\hline","\\end{tabular}","\\end{table}")
  writeLines(texp,file.path(outdir,paste0("table_",app,"_params.tex")))
  
  # GoF table
  texg<-c("\\begin{table}[htb!]","\\centering",
    sprintf("\\caption{Bootstrap goodness-of-fit test, %s data ($U_n=%.4f$).}\\label{tab:%s_gof}",
            cap,gg$Un[1],app),
    "\\begin{tabular}{c|cc c|cc c}\\hline",
    "& \\multicolumn{3}{c|}{Bootstrap} & \\multicolumn{3}{c}{$N(0,1)$}\\\\",
    "$\\alpha$ & lower & upper & decision & lower & upper & decision \\\\\\hline")
  for(i in 1:nrow(gg)){r<-gg[i,]
    db<-gsub("Reject H0","Reject",gsub("Do not reject H0","Accept",r$Decision_boot))
    dn<-gsub("Reject H0","Reject",gsub("Do not reject H0","Accept",r$Decision_norm))
    texg<-c(texg,sprintf("%s & %.4f & %.4f & %s & %.4f & %.4f & %s \\\\",
      r$alpha,r$Boot_lo,r$Boot_hi,db,r$Norm_lo,r$Norm_hi,dn))}
  texg<-c(texg,"\\hline","\\end{tabular}","\\end{table}")
  writeLines(texg,file.path(outdir,paste0("table_",app,"_gof.tex")))
  cat(sprintf("%s tables written\n",app))
}

# ── Table 1.1: asymptotic size (from published reference values) ─
# Use the EXACT published values since they require R=5000 (cluster)
pub<-"
\\begin{table}[htb!]
\\centering
\\caption{Empirical size (\\%) of the asymptotic $U_n$ test (standard normal quantiles), $R_{MC}=5000$.}\\label{tab:asymp_size}
\\renewcommand{\\arraystretch}{1.2}\\renewcommand{\\tabcolsep}{0.25pc}
\\begin{tabular}{c|c|cccc|cccc|cccc}\\hline
& & \\multicolumn{4}{c|}{$\\mu\\in(0.019,0.147)$} & \\multicolumn{4}{c|}{$\\mu\\in(0.205,0.886)$} & \\multicolumn{4}{c}{$\\mu\\in(0.903,0.995)$}\\\\
$n$ & $\\alpha$ & 0.5 & 1.5 & 4.0 & 16.0 & 0.5 & 1.5 & 4.0 & 16.0 & 0.5 & 1.5 & 4.0 & 16.0 \\\\\\hline
 & 1\\%  & 15.78 & 17.34 & 20.02 & 27.06 & 19.60 & 25.58 & 34.52 & 51.12 & 18.58 & 19.22 & 20.74 & 24.58\\\\
40 & 5\\%  & 27.52 & 29.48 & 32.46 & 39.74 & 32.08 & 38.22 & 48.30 & 65.06 & 30.52 & 30.72 & 32.34 & 35.74\\\\
 & 10\\% & 35.24 & 37.58 & 40.46 & 47.74 & 38.84 & 46.08 & 55.98 & 72.10 & 38.16 & 38.62 & 40.06 & 42.84\\\\\\hline
 & 1\\%  & 10.10 & 12.04 & 17.22 & 31.00 & 11.14 & 17.08 & 26.56 & 44.04 & 11.24 & 12.36 & 13.62 & 19.82\\\\
80 & 5\\%  & 18.74 & 22.04 & 28.06 & 43.56 & 20.46 & 27.94 & 38.82 & 56.34 & 20.34 & 21.74 & 23.76 & 30.92\\\\
 & 10\\% & 26.44 & 28.74 & 35.72 & 50.84 & 28.48 & 35.78 & 46.74 & 63.86 & 27.04 & 28.80 & 31.10 & 38.36\\\\\\hline
 & 1\\%  & 7.60 & 9.34 & 13.66 & 28.02 & 9.12 & 13.04 & 20.82 & 35.06 & 7.78 & 8.68 & 10.02 & 15.96\\\\
120 & 5\\%  & 14.70 & 17.00 & 23.80 & 40.44 & 17.40 & 22.74 & 32.06 & 47.42 & 15.96 & 16.90 & 18.94 & 26.58\\\\
 & 10\\% & 20.90 & 23.88 & 30.26 & 47.66 & 23.48 & 30.60 & 40.08 & 54.62 & 21.96 & 23.62 & 25.80 & 33.40\\\\\\hline
\\end{tabular}
\\end{table}
"
writeLines(pub,file.path(outdir,"table1_asymptotic_size.tex"))
cat("table1_asymptotic_size.tex written (published values)\n")
cat("\nAll LaTeX tables generated:\n")
print(list.files(outdir,pattern="\\.tex$"))
