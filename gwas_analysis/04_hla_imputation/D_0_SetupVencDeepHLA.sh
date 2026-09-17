#! /bin/bash
#------------------------------------------------------------------------------
#  Prepare HLA-imputation for DEEP*HLA
#------------------------------------------------------------------------------
#  ref)
#  - https://github.com/tatsuhikonaito/DEEP-HLA
#  - https://pytorch.org/get-started/previous-versions/
#
#  [Dependency of DEEP*HLA]
#  - Pytorch==1.4.0 -> refer 'Pytorch'
#       -> 1.8.1
#  - Numpy==1.17.2
#  - Pandas==0.25.1
#  - Scipy==1.3.1
#  - Argparse==1.4.0
#
#  [Pytorch==1.8.1]
# # ROCM 4.0.1 (Linux only)
# pip install torch==1.8.1+rocm4.0.1 torchvision==0.9.1+rocm4.0.1 torchaudio==0.8.1 -f https://download.pytorch.org/whl/torch_stable.html
#
# # ROCM 3.10 (Linux only)
# pip install torch==1.8.1+rocm3.10 torchvision==0.9.1+rocm3.10 torchaudio==0.8.1 -f https://download.pytorch.org/whl/torch_stable.html
#
# # CUDA 11.1
# pip install torch==1.8.1+cu111 torchvision==0.9.1+cu111 torchaudio==0.8.1 -f https://download.pytorch.org/whl/torch_stable.html
#
# # CUDA 10.2
# pip install torch==1.8.1+cu102 torchvision==0.9.1+cu102 torchaudio==0.8.1 -f https://download.pytorch.org/whl/torch_stable.html
#
# # CUDA 10.1
# pip install torch==1.8.1+cu101 torchvision==0.9.1+cu101 torchaudio==0.8.1 -f https://download.pytorch.org/whl/torch_stable.html
#
# # CPU only
# pip install torch==1.8.1+cpu torchvision==0.9.1+cpu torchaudio==0.8.1 -f https://download.pytorch.org/whl/torch_stable.html
#
#  [Pytorch==1.4.0]
# # CUDA 10.1
# pip install torch==1.4.0 torchvision==0.5.0
#
# # CUDA 9.2
# pip install torch==1.4.0+cu92 torchvision==0.5.0+cu92 -f https://download.pytorch.org/whl/torch_stable.html
#
# # CPU only
# pip install torch==1.4.0+cpu torchvision==0.5.0+cpu -f https://download.pytorch.org/whl/torch_stable.html
#
#------------------------------------------------------------------------------
set -eu

# [Main script ]
VENV_DIR=${HOME}/.venvs
VNAME=DEEP-HLA
if [ ! -d ${VENV_DIR}/${VNAME} ]; then
    wd=$(pwd)
    cd $VENV_DIR
    python3 -m venv ${VNAME}
    source ${VNAME}/bin/activate
    pip3 install \
        Numpy==1.17.2 \
        Pandas==0.25.1 \
        Scipy==1.3.1 \
        Argparse==1.4.0
    # install Pytorch==1.8.1
    pip install torch==1.8.1+cpu \
        torchvision==0.9.1+cpu \
        torchaudio==0.8.1 \
        -f https://download.pytorch.org/whl/torch_stable.html
    # install Pytorch==1.4.0
    # pip install torch==1.4.0+cpu \
    #     torchvision==0.5.0+cpu \
    #     -f https://download.pytorch.org/whl/torch_stable.html
    ${VNAME}/bin/python3 -m pip install --upgrade pip
    deactivate
    cd $wd
fi