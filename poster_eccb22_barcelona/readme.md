
~~~
convert ../out/gwas413/genome/5e-08/1000000/pltbar_x_per_variant_etissue_y_egene.py/plt.png img/pltbar_x_per_variant_etissue_y_egene.eps
convert ../out/gwas413/genome/5e-08/1000000/pltbar_x_per_variant_egene_y_etissue.py/plt.png img/pltbar_x_per_variant_egene_y_etissue.eps
convert ../out/gwas413/genome/5e-08/1000000/pltbar_x_per_variant_egene_etissue_y_gwas.py/plt.png img/pltbar_x_per_variant_egene_etissue_y_gwas.eps

convert ../out/gwas413/genome/5e-08/1000000/pltbox_x_per_rsid_y_remaptf.py/bxplt_remaptf_per_rsid_flank_0.png img/bxplt_remaptf_per_rsid_flank_0.eps

convert ../ms/fig/graphical_summary.png img/graphical_summary.eps

convert /home/gonzalez/Repositories/eqtl2gwas_pleiotropy/out/gwas413/genome/5e-08/1000000/pltbar_vep_consequence.py/missense_variant.png img/missense_variant.eps
convert /home/gonzalez/Repositories/eqtl2gwas_pleiotropy/out/gwas413/genome/5e-08/1000000/pltbar_vep_consequence.py/splice_region_variant.png img/splice_region_variant.eps
convert /home/gonzalez/Repositories/eqtl2gwas_pleiotropy/out/gwas413/genome/5e-08/1000000/pltbar_vep_consequence.py/3_prime_UTR_variant.png img/3_prime_UTR_variant.eps

latex poster_eccb22_barcelona.tex && dvipdf poster_eccb22_barcelona.dvi
~~~
