#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

hlines <- scan("hlines.list",what=numeric())
vlines <- scan("vlines.list",what=numeric())

hmax <- max(hlines)
vmax <- max(vlines)

x0 <- scan ("x0_unflagged.list",what=numeric())
x1 <- scan ("x1_unflagged.list",what=numeric())
y0 <- scan ("y0_unflagged.list",what=numeric())
y1 <- scan ("y1_unflagged.list",what=numeric())

x0_flag <- scan ("x0_flagged.list",what=numeric())
x1_flag <- scan ("x1_flagged.list",what=numeric())
y0_flag <- scan ("y0_flagged.list",what=numeric())
y1_flag <- scan ("y1_flagged.list",what=numeric())

pdf(file=args[1])
plot (0,
      xlim=c(0,vmax),
      ylim=c(0,hmax),
      xlab = "Reference",
      ylab="Query",
      main="Genome-genome alignment to identify breaks in syntenty")
abline(h=c(hlines),v=c(vlines),col="lightblue")
abline(h=c(0,hmax),v=c(0,vmax),col="green")
segments(x0,y0,x1,y1)
segments(x0_flag,y0_flag,x1_flag,y1_flag, col="red",lwd=2)
dev.off()