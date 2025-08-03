# Pancreatic Hormone Modeling

**Goal:** Building on ODE model in [Brown &amp; Tzanakakis (2023)](https://www.frontiersin.org/journals/endocrinology/articles/10.3389/fendo.2023.1212749/full) to incorporate somatostatin secretion into pancreatic islet dynamics. Simulating results in MATLAB to determine model efficacy.

**Repo Guide**: 
- Presentation.pdf: Slides used for final presentation
- Paper.pdf: Typset final report with Figures
- alphaBetaModel_perfusion_delta: Modifications of Alpha & Beta Cell ODE Model - function file for ODE model
- simulate_alphaBetaModel_perfusion_delta: Returns net signals & secretion as an array
- run_alphaBetaModel_perfusion_delta: Main file - runs simulation and produces figures

**Abstract**: It is well known that blood sugar is regulated by insulin and glucagon, two hormones produced by pancreatic islet cells. Mathematical modeling of these hormones can be a powerful tool in developing treatments for diseases such as type II diabetes (T2D). This work builds upon the model proposed by Brown and Tzanakakis (2023) by incorporating the effects of a third hormone, somatostatin, an inhibitor of both insulin and glucagon. Our model better predicted glucagon dynamics and showed that somatostatin effectively inhibits glucagon secretion in hyperglycemic conditions, with minimal impacts on insulin. This suggests that somatostatin therapy may be an effective treatment for T2D.

**Source Code**: I modified scripts from [Brown and Tzanakakis](https://github.com/aedanbrown/Paracrine-Glucose-and-Insulin-in-Glucagon-Secretion) to incorporate my new assumptions on delta-cell behavior.
