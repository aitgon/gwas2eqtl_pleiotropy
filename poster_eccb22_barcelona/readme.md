These commands run in the present folder.

~~~
mkdir -p img_tmp
convert ../out/gwas413/genome/5e-08/1000000/pltbar_x_per_variant_etissue_y_egene.py/plt.png img_tmp/pltbar_x_per_variant_etissue_y_egene.eps
convert ../out/gwas413/genome/5e-08/1000000/pltbar_x_per_variant_egene_y_etissue.py/plt.png img_tmp/pltbar_x_per_variant_egene_y_etissue.eps
convert ../out/gwas413/genome/5e-08/1000000/pltbar_x_per_variant_egene_etissue_y_gwas.py/plt.png img_tmp/pltbar_x_per_variant_egene_etissue_y_gwas.eps

convert ../out/gwas413/genome/5e-08/1000000/pltbox_x_per_rsid_y_remaptf.py/bxplt_remaptf_per_rsid_flank_0.png img_tmp/bxplt_remaptf_per_rsid_flank_0.eps

convert ../ms/fig/graphical_summary.png img_tmp/graphical_summary.eps

convert ../out/gwas413/genome/5e-08/1000000/pltbar_vep_consequence.py/missense_variant.png img_tmp/missense_variant.eps
convert ../out/gwas413/genome/5e-08/1000000/pltbar_vep_consequence.py/splice_region_variant.png img_tmp/splice_region_variant.eps
convert ../out/gwas413/genome/5e-08/1000000/pltbar_vep_consequence.py/3_prime_UTR_variant.png img_tmp/3_prime_UTR_variant.eps

convert ../out/gwas413/genome/5e-08/1000000/pltbar_davidgo.py/david_pleio_3.png img_tmp/david_pleio_3.eps
convert ../out/gwas413/genome/5e-08/1000000/pltbar_davidgo.py/david_pleio_2.png img_tmp/david_pleio_2.eps

wget https://www.univ-amu.fr/system/files/2021-02/DIRCOM-Logo-FSS.png -P img_tmp
convert img_tmp/DIRCOM-Logo-FSS.png img_tmp/DIRCOM-Logo-FSS.eps

wget https://icorc1f2.inserm.fr/content/download/82230/569063/version/1/file/logo_generique_SD.zip -P img_tmp
unzip img_tmp/logo_generique_SD.zip -d img_tmp/
convert img_tmp/logo_generique_SD/logo-generique-SD.png img_tmp/inserm.eps

wget https://tagc.univ-amu.fr/sites/tagc.univ-amu.fr/themes/webkit_theme/logo.png -O img_tmp/tagc_logo.png
convert img_tmp/tagc_logo.png img_tmp/tagc_logo.eps

latex poster_eccb22_barcelona.tex && dvipdf poster_eccb22_barcelona.dvi
~~~
