while read line1
do
   read line2
   PYTHONUNBUFFERED=1 python fold_SM.py -f $line1 -s $line2 -k 3 -a 3 -e -1 -u False -d 2 2>/dev/null
#    echo "-f $line1 -s $line2"
# "PYTHONUNBUFFERED=1 ~/tools/miniconda3/envs/py2/bin/python fold_SM.py {} -k 3 -a 3 -e -1 -u False -d 2 2>/dev/null"
done 
# | parallel --keep-order "module load python/python-2.7.11 &&  PYTHONUNBUFFERED=1 python2 fold_SM.py {} -k 3 -a 3 -e -1 -u False -d 2 2>/dev/null"