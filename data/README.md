# Input data contract

Provide smile_period_means.npz with no object or pickle arrays. Required arrays:

- model_names: 12 distinct model names, stored as a Unicode string array.
- latitude and longitude: 1D coordinates for a regular latitude/longitude grid.
- mask: boolean study-region mask of shape (latitude, longitude).
- rx_his, rx_future: Rx1day climatologies in mm/day.
- q_his, q_future: validated 2 m saturation specific humidity climatologies in kg/kg.
- tas_his, tas_future: temperature climatologies with consistent units (K or Celsius).

The six fields have shape (12, latitude, longitude) and identical coordinate alignment. Each model must contribute one mean of its verified SMILE members. Historical period is 1985-2014; future period is 2070-2099 under SSP5-8.5/RCP8.5. Record model/member counts and IDs, scenario, grids, masks, units, averaging order, source checksums, and preprocessing outside this array contract.

Optionally provide observed_rx1day_climatology.npz with rx1day_climatology in mm/day, latitude, and longitude. It must be the verified 1985-2014 observation field, already on the identical model grid. The runner does not merge observation products or select a regridding method.

No data are supplied in this package. Missing primary inputs cause an error; no substitute or random data are generated.
