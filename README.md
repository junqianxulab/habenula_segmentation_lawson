# Improved implementation<sup>[1](#kim)</sup> of Lawson's geometric habenula segmentation<sup>[2](#lawson)</sup>

---

### Depencency

##### with Docker
* Docker
* X11 server
    * linux: Most desktop linux distributions contain X11 server. If not, install Xorg and/or openbox.
    * macOS: XQuartz
    * Windows: VcXsrv, Xming, MobaXterm, ... but it has not been tested.

##### without Docker
* python 2.7
* numpy
* nibabel
* pyside

---

### Input
* T1w_acpc_dc.nii.gz (AC-PC-aligned T1-weighted image)

---

### Installation and usage

##### with Docker
* Make sure you are in docker group. If not, run `sudo usermod -aG docker $USER` and logout/login.
* macOS users should check "Allow connection from network clients" option in the Security tab in Xquartz Preferences.

* Build docker: in `src` directory, run `docker build -t hb_lawson .`
* `src/run_hb_lawson_docker.sh seg T1w_acpc_dc.nii.gz my_habenula_segmentation_name` : open segmentation program
* `src/run_hb_lawson_docker.sh to_nifti T1w_acpc_dc.nii.gz my_habenula_segmentation_name_lawson_hb.csv` : convert csv file to nifti file

* It is convenient to copy `src/run_hb_lawson_docker.sh` file to one of your PATH environment directory.
* macOS user note: if detected IP address is wrong, modify `src/run_hb_lawson_docker.sh` file:line 10.
* Windows user note: this has not been tested on Windows.

##### without Docker
* `python hb_segmentation_lawson.py T1w_acpc_dc.nii.gz my_habenula_segmentation_name` : open segmentation program
* `python hb_lawson_to_nifti.py T1w_acpc_dc.nii.gz my_habenula_segmentation_name_lawson_hb.csv` : convert csv file to nifti file

---

### Reference
<a name="Kim">[1]</a> Kim et al., Reproducibility of myelin content‐based human habenula segmentation at 3 Tesla. Human Brain Mapping. 2018; 39: 3058-3071 https://doi.org/10.1002/hbm.24060

<a name="lawson">[2]</a> Lawson et al., Defining the habenula in human neuroimaging studies. Neuroimage. 2013; 64: 722–727 https://doi.org/10.1016/j.neuroimage.2012.08.076

