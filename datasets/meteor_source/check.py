from obspy import read
st = read("path_to/4_sonic_koblenz_waveforms_response_removed.mseed")
st_filt = st.copy().filter("bandpass",corners=4,freqmin=1,freqmax=10,zerophase=False)
st_filt.select(station="DEP52").plot(outfile='path_to/check.png')
