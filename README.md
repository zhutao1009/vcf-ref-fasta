# vcf-ref-fasta
vcf2fa1.py and vcf2fa2.py are 2 script with same function.  
vcf2fa1.py is less requirement(numpy), but can only conver a vcf file with small reference genome, 100M maybe.  
vcf2fa2.py depend on numpy and [pyfasta](https://github.com/brentp/pyfasta.git), which is more efficient can processing big reference genome, over 1G.  


usage: vcf2fa1.py [-h] -v VCF -r REF [-k KEEP]  
usage: vcf2fa2.py [-h] -v VCF -r REF [-k KEEP]
-h help message  
-v your vcf file  
-r your reference squence  

zhutao@cau.edu.cn for help
