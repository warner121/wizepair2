#!/bin/sh
# Copyright 2022 Google LLC
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

#sudo apt-get install software-properties-common
#sudo add-apt-repository ppa:deadsnakes/ppa
sudo apt-get update
sudo apt-get -y install python3
sudo apt-get -y install python3-distutils
sudo apt-get -y install wget

export PATH="$PATH:/local/bin"
wget -nv https://bootstrap.pypa.io/get-pip.py
python3 get-pip.py

dir=/mnt/share
infile=$dir/mmp_input/mmp-`printf %0.12d $BATCH_TASK_INDEX`.json.gz
outfile=$dir/mmp_output/mmp-`printf %0.12d $BATCH_TASK_INDEX`.json.gz

mkdir -p $dir/mmp_output

rm -rf wizepair2
git clone https://github.com/warner121/wizepair2.git 
cd wizepair2

pip3 install --root-user-action=ignore -r requirements.txt
#python3 unit_tests.py

python3 batch/mmp_batch.py $infile $outfile 
