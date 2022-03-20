
run.lw.estimation <- function(log.output, inflation, real.interest.rate, nominal.interest.rate,
                               a3.constraint = NA, b2.constraint = NA, run.se = TRUE) {
   
    out.stage1 <- rstar.stage1(log.output,
                               inflation,
                               b2.constraint)

    lambda.g <- median.unbiased.estimator.stage1(out.stage1$potential.smoothed)

    out.stage2 <- rstar.stage2(log.output,
                               inflation,
                               real.interest.rate,
                               lambda.g,
                               a3.constraint,
                               b2.constraint)

    lambda.z <- median.unbiased.estimator.stage2(out.stage2$y, out.stage2$x)

    out.stage3 <- rstar.stage3(log.output,
                               inflation,
                               real.interest.rate,
                               nominal.interest.rate,
                               lambda.g,
                               lambda.z,
                               a3.constraint,
                               b2.constraint,
                               run.se)

    return(list(out.stage1=out.stage1,out.stage2=out.stage2,out.stage3=out.stage3,
                lambda.g=lambda.g,lambda.z=lambda.z))
}
