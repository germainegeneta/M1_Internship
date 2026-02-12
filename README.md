# M1_Internship
R scripts for cross-correlation mapping of mosquito population responses and climate variables.

1. Transforming to DwC format - harmonizes entomological datasets to comply with DarwinCore (DwC) standards. Outputs are data tables following DwC sampling-event, occurrence, and metadata dataset requirements.

2. OpenMeteo data retrieval - downloads historical climate data from OpenMeteo via r/openmeteo. Outputs are data tables containing climate measures for the specified time and location.

3. Cross-correlation mapping using climwin - generates CCMs via r/climwin. Outputs are site-scale and regional-scale maps computed using beta coefficients and deltaAICc. 