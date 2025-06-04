# Quantifying Supply Chain Drivers of Epidemic Aedes-borne Virus Risk for Systemic Prevention

## Code used to inform study design
This repository contains the code used to inform the study design for the entomological and epidemiological aspects of our UKRI application using power analysis by simulation.

## Contents

- **MosquitoAbundanceSampling.R**: power analysis for study design to detect differences in Aedes abundance between land use types.

- **Seroprevalence.R**: power analysis for cross-sectional surveys to estimate differences in seroprevalence between high and low risk areas.

  - **ForceOfInfection.R**: power analysis to inform study desing to estimate time-varying force of infection.

You will need to install the GLMMmisc which can be found at devtools::install_github("pcdjohnson/GLMMmisc")
