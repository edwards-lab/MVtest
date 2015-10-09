
import numpy
import scipy.stats
from pygwas.pheno_covar import PhenoCovar
import pygwas.standardizer
import sys

class Standardizer(pygwas.standardizer.StandardizedVariable):
    """Optional plugin object that can be used to standardize covariate and
       phenotype data.

       Many algorithms require that input be standardized in some way in order
       to work properly, however, rescaling the results is algorithm specific.
       In order to facilitate this situation, application authors can
       write up application specific Standardization objects for use with
       the data parsers.

       """

    def __init__(self, pc):
        super(Standardizer, self).__init__(pc)


    def standardize(self):
        """Standardize the variables within a range [-1.0 and 1.0]

        This replaces the local copies of this data. When it's time to
        scale back, use destandardize from the datasource for that.

        """


        nonmissing = numpy.invert(self.missing)
        nmcount = numpy.sum(nonmissing)

        self.covariates = []
        self.phenotypes = []

        for idx in range(0, self.covar_count):
            x = self.datasource.covariate_data[idx]
            nonmissing = x != PhenoCovar.missing_encoding
            mx = numpy.mean(x[nonmissing])
            sx = numpy.std(x[nonmissing])
            self.covariates.append((x-mx)/sx)

        #for idx in range(0, self.pheno_count):
        y = self.datasource.phenotype_data[self.idx]
        nonmissing = y != PhenoCovar.missing_encoding
        my = numpy.mean(y[nonmissing])
        sy = numpy.std(y[nonmissing])
        self.phenotypes.append((y-my)/sy)

    def destandardize(self, estimates, se, **kwargs):
        """Revert the betas and variance components back to the original scale.

        """
        pvalues = kwargs["pvalues"]
        v = kwargs["v"]
        nonmissing=kwargs["nonmissing"]
        pheno = self.datasource.phenotype_data[self.idx][self.datasource.phenotype_data[self.idx] != PhenoCovar.missing_encoding]
        covariates = []
        mmx = []
        ssx = []
        a = [1,0]

        for c in self.datasource.covariate_data:
            covariates.append(c[c!=PhenoCovar.missing_encoding])
            mmx.append(numpy.mean(covariates[-1]))
            ssx.append(numpy.std(covariates[-1]))
            a.append(-mmx[-1]/ssx[-1])
        if len(mmx) < 1:
            mmx = 1
            ssx = 1
        else:
            mmx = numpy.array(mmx)
            ssx = numpy.array(ssx)
        covariates = numpy.array(covariates)
        ssy = numpy.std(pheno)
        mmy = numpy.mean(pheno)
        # Quick scratch pad for Chun's new destandardization
        ccount = len(covariates) + 2
        a=numpy.array(a)

        meanpar = list([mmy + ssy * (estimates[0] - numpy.sum(estimates[2:ccount]*mmx/ssx)),
                    estimates[1] * ssy])


        meanse = [ssy*numpy.sqrt(a.dot(v[0:ccount, 0:ccount]).dot(a.transpose())),
                  ssy * se[1]]
        varpar = [2*numpy.log(ssy) + estimates[ccount] - numpy.sum(estimates[ccount+2:]*mmx/ssx),
                  estimates[ccount+1]]
        varse = [numpy.sqrt(a.dot(v[ccount:, ccount:]).dot(a.transpose())),
                 se[ccount+1]]
        if ccount > 2:
            meanpar += list(estimates[2:ccount] * ssy / ssx)
            meanse += list(ssy * se[2:ccount]/ssx)
            varpar += list(estimates[ccount+2:]/ssx)
            varse += list(se[ccount+2:]/ssx)
        pvals = 2*scipy.stats.norm.cdf(-numpy.absolute(numpy.array(meanpar+varpar)/numpy.array(meanse+varse)))



        return meanpar+varpar, meanse + varse, pvals