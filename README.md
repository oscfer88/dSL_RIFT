# Statistical Learning of Distractor Suppression Down-regulates Pre-Stimulus Neural Excitability in Early Visual Cortex
This repository contains the experimental code used to run the experiment, all the scripts used to analyse the data, and the questionnaire used to assess participants' awareness of the statistical learning manipulation.
<br /><br />

### **The study is pubblished on Jornal of Neuroscience. Link [here!](https://doi.org/10.1523/JNEUROSCI.1703-22.2022)**
<br />

## Abstract:
Visual attention is highly influenced by past experiences. Recent behavioral research has shown that expectations about the spatial location of distractors within a search array are implicitly learned, with expected distractors becoming less interfering. Little is known about the neural mechanism supporting this form of statistical learning. Here we used magnetoencephalography (MEG) to measure human brain activity to test whether proactive mechanisms are involved in the statistical learning of distractor locations. Specifically, we used a new technique called rapid invisible frequency tagging (RIFT) to assess neural excitability in early visual cortex during statistical learning of distractor suppression, while concurrently investigating the modulation of posterior alpha-band activity (8-12 Hz). Male and female human participants performed a visual search task in which a target was occasionally presented alongside a color-singleton distractor. Unbeknown to the participants, the distracting stimuli were presented with different probabilities across the two hemifields. RIFT analysis showed that early visual cortex exhibited reduced neural excitability in the pre-stimulus interval at retinotopic locations associated with higher distractor probabilities. In contrast, we did not find any evidence of expectation-driven distractor suppression in alpha-band activity. These findings indicate that proactive mechanisms of attention are involved in predictive distractor suppression and that these mechanisms are associated with altered neural excitability in early visual cortex. Moreover, our findings indicate that RIFT and alpha-band activity might subtend different and possibly independent attentional mechanisms.
<br />

## Experimental paradigm:
<img src="https://github.com/oscfer88/dSL_RIFT/blob/main/_pics/01%20paradigm.png?raw=true" width=50% height=50%>
<br />
The visual search task. Each trial started with a fixation dot. Afterwards, a placeholder screen with four Gabor patches was displayed. Then the search array was presented, and participants had to report whether the target was tilted to the left or right. In 66% of the trial, the target was presented together with a color-singleton distractor. During the whole RIFT periods, the two stimuli at the bottom were showed flickering at a broadband signal (RIFT). Statistical learning manipulation. The distractor was presented more often in one hemifield (75%) than the other (25%). To counterbalance the distractor probability conditions across the two hemifields, each participant took part in two experimental sessions (indicated as Day 1 and Day 2), with the statistical learning manipulation swapped after the first session.
<br /><br />

## Rapid Invisible Frequency Tagging
<img src="https://github.com/oscfer88/dSL_RIFT/blob/main/_pics/02%20freq_tagging.png?raw=true" width=80% height=80%>
<br />
During the frequency tagging period, the two stimuli at the bottom of the screen (indicated by the dotted orange and blue circles) were flickered at very high frequencies (> 50 Hz) to generate a steady-state response in visual cortex. Examples of two uncorrelated (i.e., left and right stimulation signals) broadband tagging signals are shown in the figure above, as well as the power spectra of the broadband tagging signals used in this study.
<br /><br />

## Software and Toolboxes:
- MATLAB
- [FieldTrip](https://www.fieldtriptoolbox.org/)
- [computeCohen_d](https://uk.mathworks.com/matlabcentral/fileexchange/62957-computecohen_d-x1-x2-varargin)
- [bayesFactor](https://github.com/klabhub/bayesFactor)
- [violinplot](https://github.com/bastibe/Violinplot-Matlab)
