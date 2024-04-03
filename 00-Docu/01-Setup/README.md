This document will cover the installation of the necessary dependencies for the project

1. Install python. The version i tested on is 3.10.12
2. Create .venv (highly recommended)
3. Run Makefile:
    On new terminal, make sure that it is preceded by (.venv) if used
    
    make

4. Install cugraph. For more information reffer to (https://docs.rapids.ai/install)
    In my case, the command used was:
    
    pip install \
        --extra-index-url=https://pypi.nvidia.com \
        cudf-cu12==24.2.* dask-cudf-cu12==24.2.* cuml-cu12==24.2.* \
        cugraph-cu12==24.2.* cuspatial-cu12==24.2.* cuproj-cu12==24.2.* \
        cuxfilter-cu12==24.2.* cucim-cu12==24.2.* pylibraft-cu12==24.2.* \
        raft-dask-cu12==24.2.*



For the future. Makefile would really benefit from specifying dependencie versions, as i havent done it yet, i've included the installation trace (in order to check each of the versions)


------------------------------------------------------------------------------------

make
Requirement already satisfied: matplotlib in ./.venv/lib/python3.10/site-packages (from -r requirements.txt (line 6)) (3.8.3)
Requirement already satisfied: seaborn in ./.venv/lib/python3.10/site-packages (from -r requirements.txt (line 7)) (0.13.2)
Requirement already satisfied: numpy in ./.venv/lib/python3.10/site-packages (from -r requirements.txt (line 8)) (1.24.4)
Requirement already satisfied: faiss-cpu in ./.venv/lib/python3.10/site-packages (from -r requirements.txt (line 11)) (1.8.0)
Collecting faiss-gpu (from -r requirements.txt (line 12))
  Using cached faiss_gpu-1.7.2-cp310-cp310-manylinux_2_17_x86_64.manylinux2014_x86_64.whl.metadata (1.4 kB)
Collecting anndata (from -r requirements.txt (line 13))
  Using cached anndata-0.10.5.post1-py3-none-any.whl.metadata (6.5 kB)
Requirement already satisfied: contourpy>=1.0.1 in ./.venv/lib/python3.10/site-packages (from matplotlib->-r requirements.txt (line 6)) (1.2.0)
Requirement already satisfied: cycler>=0.10 in ./.venv/lib/python3.10/site-packages (from matplotlib->-r requirements.txt (line 6)) (0.12.1)
Requirement already satisfied: fonttools>=4.22.0 in ./.venv/lib/python3.10/site-packages (from matplotlib->-r requirements.txt (line 6)) (4.49.0)
Requirement already satisfied: kiwisolver>=1.3.1 in ./.venv/lib/python3.10/site-packages (from matplotlib->-r requirements.txt (line 6)) (1.4.5)
Requirement already satisfied: packaging>=20.0 in ./.venv/lib/python3.10/site-packages (from matplotlib->-r requirements.txt (line 6)) (24.0)
Requirement already satisfied: pillow>=8 in ./.venv/lib/python3.10/site-packages (from matplotlib->-r requirements.txt (line 6)) (10.2.0)
Requirement already satisfied: pyparsing>=2.3.1 in ./.venv/lib/python3.10/site-packages (from matplotlib->-r requirements.txt (line 6)) (3.1.2)
Requirement already satisfied: python-dateutil>=2.7 in ./.venv/lib/python3.10/site-packages (from matplotlib->-r requirements.txt (line 6)) (2.9.0.post0)
Requirement already satisfied: pandas>=1.2 in ./.venv/lib/python3.10/site-packages (from seaborn->-r requirements.txt (line 7)) (1.5.3)
Collecting array-api-compat (from anndata->-r requirements.txt (line 13))
  Using cached array_api_compat-1.5-py3-none-any.whl.metadata (16 kB)
Requirement already satisfied: exceptiongroup in ./.venv/lib/python3.10/site-packages (from anndata->-r requirements.txt (line 13)) (1.2.0)
Collecting h5py>=3 (from anndata->-r requirements.txt (line 13))
  Using cached h5py-3.10.0-cp310-cp310-manylinux_2_17_x86_64.manylinux2014_x86_64.whl.metadata (2.5 kB)
Collecting natsort (from anndata->-r requirements.txt (line 13))
  Using cached natsort-8.4.0-py3-none-any.whl.metadata (21 kB)
Requirement already satisfied: scipy>1.4 in ./.venv/lib/python3.10/site-packages (from anndata->-r requirements.txt (line 13)) (1.12.0)
Requirement already satisfied: pytz>=2020.1 in ./.venv/lib/python3.10/site-packages (from pandas>=1.2->seaborn->-r requirements.txt (line 7)) (2024.1)
Requirement already satisfied: six>=1.5 in ./.venv/lib/python3.10/site-packages (from python-dateutil>=2.7->matplotlib->-r requirements.txt (line 6)) (1.16.0)
Using cached faiss_gpu-1.7.2-cp310-cp310-manylinux_2_17_x86_64.manylinux2014_x86_64.whl (85.5 MB)
Using cached anndata-0.10.5.post1-py3-none-any.whl (121 kB)
Using cached h5py-3.10.0-cp310-cp310-manylinux_2_17_x86_64.manylinux2014_x86_64.whl (4.8 MB)
Using cached array_api_compat-1.5-py3-none-any.whl (39 kB)
Using cached natsort-8.4.0-py3-none-any.whl (38 kB)
Installing collected packages: faiss-gpu, natsort, h5py, array-api-compat, anndata
Successfully installed anndata-0.10.5.post1 array-api-compat-1.5 faiss-gpu-1.7.2 h5py-3.10.0 natsort-8.4.0
