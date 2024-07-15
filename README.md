# Stata-code
The .ado files are Stata programs, the .hlp files are Stata help files. Some of the programs require that you have other programs installed. If you encounter an "unrecognized program: <program name>" error in STATA, you can type "findit <program name>" to try to find the missing program on line. You can install any of the following programs individually by copying the relevant .ado and .hlp files to your ado directory (see here). You can find the directory by typing "personal" into the Stata command line.

See the UCLA Academic Technology Services' excellent statistical computing website http://www.ats.ucla.edu/stat/overview.htm for a variety of examples and references. Also see the programs on the Michigan Accounting PhD web site.

<table xmlns="http://www.w3.org/1999/xhtml" border="1" width="757" data-table-local-id="table-3">
  <caption>
    <div align="left">
      <b>Estimating equations:</b>
    </div>
  </caption>
  <colgroup>
    <col width="140" align="left">
      <col width="360" align="left">
  </colgroup>
  <thead style="border-bottom:1px solid black;">
    <tr>
      <th scope="col">File(s)</th>
      <th scope="col">Description</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td align="left" valign="top">
        bking2.ado
      </td>
      <td align="left" valign="top">Baxter-King (1995) bandpass filter for time series.</td>
    </tr>
    <tr>
      <td align="left" valign="top">
        xtfmbJ.ado, xtfmbJ.hlp
      </td>
      <td align="left" valign="top">Fama-MacBeth (1973) regressions with options to weight by number of observations as proxy for precision of the years' estimates and an option to use a Newey-West correction for serial correlation in coefficient estimates.</td>
    </tr>
    <tr>
      <td align="left" valign="top">
        cgmreg.ado, cgmreg.hlp, cgmregF.ado, cgmregF.hlp, cgmlogit.ado, cgmlogit.hlp, cgmflogit.ado, cgmflogit.hlp, cgmwildboot.ado, cgmwildboot.hlp
      </td>
      <td align="left" valign="top">
        Implementation of various estimation commands with multi-way clustered standard errors as in Cameron, Gelbach and Miller (2010 <i>Journal of Business and Economic Statistics</i>). Also see Petersen (2009 <i>Review of Financial Studies</i>). Commands include linear regression (cgmreg), linear regression with Fieller (1954) confidence intervals on coefficient ratios (cgmregF), logit (cgmlogit), fractional logit (cgmflogit) and regression with bootstrapped p-values (cgmwildboot - See Cameron, Gelbach and Miller, 2008 <i>Review of Economics and Statistics</i>).
        <p>
          - cgmwildboot updated 2013-09-06
          <br />
          - Updated 2014-02-19 to deal with dropped variables (e.g., fixed effects dropped due to collinearity)
          <br />
          - Updated cgmwildboot 2014-03-26 to fix error introduced by prior fix for dropped variables
          <br />
          - Updated cgmreg 2014-07-03 to disallow pweight
          <br />
          - Updated cgmwildboot 2015-03-10 to fix p-values (e.g., program looks at upper-tail for p-values of positive coefficients, and this was giving p-values greater than one when the coefficient was in the bottom half of the bootstrap distribution)
          <br />
        </p>
      </td>
    </tr>
    <tr>
      <td align="left" valign="top">
       regFieller.ado, regFieller.hlp
      </td>
      <td align="left" valign="top">
        OLS regressions with confidence intervals for ratios of regression coefficients based on Fieller's theorem (more robust than delta method)
        <br />
        - Updated 2014-02-19 to deal with dropped variables (e.g., fixed effects dropped due to collinearity)
        <br />
      </td>
    </tr>
    <tr>
      <td align="left" height="44" valign="top">
       shumhaz.ado, shumhaz.hlp
      </td>
      <td align="left" valign="top">
        Shumway (2001) hazard model estimates, which uses a standard logit routine and corrects the chi-squared statistics for the average number of observations per cross-sectional unit.
        <br />
        <br />
        - Updated 2014-02-10 to include option
        <span style="font-family:courier new,monospace">lroc</span>
        to report area under ROC curve
        <br />
      </td>
    </tr>
    <tr>
      <td align="left" valign="top" width="145">
        <a>mishkin.ado</a>
      </td>
      <td align="left" valign="top" width="589">Implementation of Mishkin (1983) rational expectations tests. See comments in the ado file for syntax.</td>
    </tr>
    <tr>
      <td align="left" valign="top">
       vuong.ado, vuong.hlp
      </td>
      <td align="left" valign="top">Computes Vuong (1989 Econometrica) test of two non-nested regressions as implemented and described in Dechow (1994 Journal of Accounting and Economics). Syntax is "vuong mod1 mod2" where mod1 and mod2 are stored regression results.</td>
    </tr>
  </tbody>
</table>


<table xmlns="http://www.w3.org/1999/xhtml" border="1" width="757" data-table-local-id="table-5">
  <caption>
    <div align="left"><b>Utilities:</b></div>
  </caption>
  <colgroup>
    <col width="140" align="left">
      <col width="360" align="left">
  </colgroup>
  <thead style="border-bottom:1px solid black;">
    <tr>
      <th scope="col">File(s)</th>
      <th scope="col">Description</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td align="left" valign="top" width="145">
        <a>bpagan2.ado</a>
      </td>
      <td align="left" valign="top" width="589">Breusch-Pagan test of heteroskedasticity.</td>
    </tr>
    <tr>
      <td align="left" valign="top">
        <a>collapseJ.ado</a>
      </td>
      <td align="left" valign="top">Variant of STATA's collapse command that preserves variable labels.</td>
    </tr>
    <tr>
      <td align="left" valign="top">
        ffind.ado, ffind.hlp
      </td>
      <td align="left" valign="top">Creates Fama-French industry classifications based on SIC codes</td>
    </tr>
    <tr>
      <td align="left" valign="top">
        <a>prchgJ.ado</a>
      </td>
      <td align="left" valign="top">Computes estimated change in probabilities for user-specified changes in variables following logit/probit estimation.</td>
    </tr>
    <tr>
      <td align="left" valign="top">
        <a>prvalueJ.ado</a>
      </td>
      <td align="left" valign="top">Computes predicted probabilties and confidence intervals following logit/probit estimation.</td>
    </tr>
    <tr>
      <td align="left" height="42" valign="top">
        <a>truncateJ.ado</a>
      </td>
      <td align="left" valign="top">Creates an indicator variable that                  identifies observations for which a set of                  variables are within specified cutoffs.</td>
    </tr>
    <tr>
      <td align="left" valign="top">
        <a>winsorizeJ.ado</a>
      </td>
      <td align="left" valign="top">Winsorizes variables, allowing, for example, within year truncation points</td>
    </tr>
    <tr>
      <td align="left" height="24" valign="top">
        <a>xtileJ.ado</a>
      </td>
      <td align="left" valign="top">Computes a variable that contains quantiles. The command is similar to STATA's xtile command except that xtileJ allows 'by' groups. This program runs much more quickly than the xtile2 command that allows 'by' groups.</td>
    </tr>
  </tbody>
</table>

<table xmlns="http://www.w3.org/1999/xhtml" border="1" width="757" data-table-local-id="table-4">
  <caption>
    <div align="left">
      <b>Descriptive statistics:</b>
    </div>
  </caption>
  <colgroup>
    <col width="140" align="left">
      <col width="360" align="left">
  </colgroup>
  <thead style="border-bottom:1px solid black;">
    <tr>
      <th scope="col">File(s)</th>
      <th scope="col">Description</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td align="left" valign="top">
        <a>corrtbl.ado</a>
      </td>
      <td align="left" valign="top">Creates correlation matrix for a list of variables with Pearson correlations below the diagonal and Spearman correlations above.</td>
    </tr>
    <tr>
      <td align="left" valign="top">
        <a>counttex.ado</a>
      </td>
      <td align="left" valign="top">Generates LaTeX code for a frequency table.</td>
    </tr>
    <tr>
      <td align="left" valign="top">
        <a>summstattable.ado</a>
      </td>
      <td align="left" valign="top">Generates LaTeX code for a table of summary statistics for a set of variables grouped by a binary indicator. Also reports t-tests for difference in means and chi-squared test for difference in medians.</td>
    </tr>
    <tr>
      <td align="left" valign="top" width="145">
        <a>ttesttable.ado</a>
      </td>
      <td align="left" valign="top" width="589">Generates LaTeX code for a table that displays means and t-tests for a set of variables grouped by a binary indicator.</td>
    </tr>
  </tbody>
</table>
