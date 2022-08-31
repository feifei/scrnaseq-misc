#!/bin/bash

PATH="/Users/feifei/miniconda3/bin:/Users/feifei/miniconda3/condabin:/Users/feifei/Tools/pyvenn:/Users/feifei/Tools/python_lib:/usr/local/sbin:/usr/local/bin:/sbin:/usr/sbin:/bin:/usr/bin"

source /Users/feifei/.bash_profile

conda activate cellxgene

cellxgene_dir="/Users/feifei/Projects/cellxgene"
cd $cellxgene_dir
cellxgene launch --host 0.0.0.0 -p 5006 --annotations-file cellxgene_anno/y.csv h5ad/macro_ss2.cellxgene.h5ad &
cellxgene launch --host 0.0.0.0 -p 5007 --annotations-file cellxgene_anno/x.csv h5ad/tabula_senis_facs.cellxgene.h5ad &
cellxgene launch --host 0.0.0.0 -p 5008 --annotations-file cellxgene_anno/y.csv h5ad/oprescu.cellxgene.h5ad &
cellxgene launch --host 0.0.0.0 -p 5009 --annotations-file cellxgene_anno/z.csv h5ad/facs_large_intestine.cellxgene.h5ad &
cellxgene launch --host 0.0.0.0 -p 5010 --annotations-file cellxgene_anno/x.csv h5ad/tabula_senis_10x.cellxgene.h5ad &
cellxgene launch --host 0.0.0.0 -p 5011 --annotations-file cellxgene_anno/x.csv h5ad/macro_ss3.cellxgene.h5ad &
cellxgene launch --host 0.0.0.0 -p 5012 --annotations-file cellxgene_anno/x.csv h5ad/macro_int.cellxgene.h5ad &

