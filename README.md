<article class="markdown-body entry-content" itemprop="mainContentOfPage"><h1><a id="user-content-myspm" class="anchor" href="#myspm" aria-hidden="true"><span class="octicon octicon-link"></span></a>myspm</h1>

<p>SPM (Statistical Parametric Mapping) is a well-known MATLAB toolbox for neuroimaging analysis. Most of analyses are programable using the MATLAB batch editor including generating result reports with tables (“Results Report”).</p>

<p>But plotting orthogonal sections of a structural image, with pretty suprathreshold blobs, overlaid is not available from the MATLAB batch editor. When you have many models and contrast to test, exploring all different models and contrasts require A TONS of mouse clicking, which might cause carpal tunnel syndrome.</p>

<p>In order to protect carpal tunnels of neuroimaging researchers with a lot of models (or just for lazy guys like me), it is possible to create a script that generates orthogonal sections at each local maximum using ‘xSPM’, ‘hReg’, and ‘TabDat’ variables (which are generated when Results Report is run), and ‘spm_sections.m’ and ‘spm_orthviews.m’ functions.</p>

<ul>
<li><code>myspm_glm.m</code> performs statistical inference using SPM's GLM modules.</li>
<li><code>myspm_results.m</code> generates results reports of maximal intensity projection(MIP or 'glass brain') and orthogonal slices view and a scatter plot for each suprathreshold cluster.</li>
<li><code>myspm_graph.m</code> creates a scatter plot and a linear regression line.</li>
<li><code>myfsl_atlasquery.m</code> returns the probablistic atals information of a given MNI coordiante using Harvard-Oxford cortical/subcortical atlases, cerebellar atlas, and JHU-ICBM DTI atlas included in FSL ($FSLDIR required to be valid).</li>
</ul>
</article>
  </div>

</div>

