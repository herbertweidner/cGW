This is a collection of MATLAB programs for searching for phase-modulated signals in any long-term records. The records can originate from sensors such as barometers or gravimeters and should span periods of at least several months. Whether and how data gaps or jumps are eliminated is up to each user.

The programs search for continuous signals of (almost) constant frequency. Signals from sources outside the solar system are always phase-modulated. This means that the reception frequency is greater than the mean value (blueshift) for six months and then smaller (redshift). It is often possible to measure phase modulations with other time durations and unknown causes.

The technical background is described in ‘The Search for Gravitational Waves: Fundamentals of Reception Technology’ (https://vixra.org/abs/2311.0020).

The file yDWD.mat contains air pressure measurements from approximately 50 DWD weather stations over a period of 20 years (January 1, 2000, to December 31, 2019). One measurement per hour. The raw data are added without correction; data gaps can lead to jumps. These are insignificant when searching for cGW with periods of several days.

