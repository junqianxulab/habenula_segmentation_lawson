FROM python:2.7

RUN apt-get -y update
RUN apt-get install -y libqt4-dev cmake xvfb
RUN pip install numpy nibabel pyside

COPY hb_lawson_to_nifti.py /usr/local/bin/
COPY hb_segmentation_lawson.py /usr/local/bin/
RUN chmod ugo+rx /usr/local/bin/hb_lawson_to_nifti.py /usr/local/bin/hb_segmentation_lawson.py
ENV QT_X11_NO_MITSHM=1

VOLUME ["/data"]
WORKDIR /data
CMD bash

