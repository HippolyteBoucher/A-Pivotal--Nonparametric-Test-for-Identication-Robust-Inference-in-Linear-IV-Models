# A Pivotal Nonparametric Test for Identication-Robust Inference in Linear IV Models


This repository contains the latest version of the working paper titled ''A Pivotal Nonparametric Test for Identication Robust Inference in Linear IV Models'' by Hippolyte Boucher. Replication files for both the simulations and application are at disposal. This working paper follows the CC-BY-NC-SA-4.0 license, it can be distributed for noncommercial purposes only, and only so long as attribution is given to the creator. Contact the writer at **Hippolyte.Boucher@outlook.com**.


## Abstract 

In linear models with endogenous regressors it is well-known that weak instruments (IVs) bias the 2 Stage Least Squares (2SLS) and other k-class IV estimators and make standard Gaussian confidence intervals invalid. Inference can still be performed by inverting tests, however there are no known method to account for a non-linear first stage except Antoine and Lavergne (2019). Their method requires simulations of the distribution of the test statistic under the null which makes it difficult to apply when sample size is moderate to large. For the above reasons I build a pivotal test statistic based on a score of integrated conditional moments which allows to easily infer on the model's structural parameters regardless of instruments' strength and the shape of the first stage conditional mean. For heteroskedastic or independent and identically distribution data with normal or non-normal errors I prove that the test is valid regardless of the degree of identification of the structural parameter of interest, and also prove that the test is consistent as long if the parameter of interest is at least semi-strongly identified. I compare the performances of the test against competing ones and revisit the effect of education on wage using Angrist and Krueger (1991) data and prove that it is strictly positive.

Keywords: Weak Instruments, Hypothesis Testing,  Semiparametric Model

JEL Codes: C12, C13, C14

## References

Bertille Antoine and Pascal Lavergne (2019), ["Identification-Robust Nonparametric Inference in a Linear IV Model"][1], *Simon Fraser University Discussion Papers*

Joshua D. Angrist and  Alan B. Krueger (1991), ["Does Compulsory School Attendance Affect Schooling and Earnings?"][2], *The Quarterly Journal of Economics*, 106 (4), 979–1014

[1]: https://ideas.repec.org/p/sfu/sfudps/dp19-02.html
[2]: https://www.jstor.org/stable/2937954
