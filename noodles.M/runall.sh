for i in N*; do echo $i; cd $i; Rscript set.and.run.test.R; cd .. ;done
