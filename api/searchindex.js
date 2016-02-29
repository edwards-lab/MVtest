Search.setIndex({envversion:46,filenames:["index","modules","modules/meanvar","modules/pygwas"],objects:{"":{meanvar:[2,0,0,"-"],pygwas:[3,0,0,"-"]},"meanvar.mv_esteq":{MeanVarEstEQ:[2,4,1,""],RunAnalysis:[2,4,1,""],RunMeanVar:[2,4,1,""]},"meanvar.mvresult":{MVResult:[2,3,1,""]},"meanvar.mvresult.MVResult":{beta_pvalues:[2,1,1,""],beta_stderr:[2,1,1,""],betas:[2,1,1,""],chr:[2,1,1,""],covar_labels:[2,1,1,""],eff_alcount:[2,1,1,""],lmpv:[2,1,1,""],maf:[2,1,1,""],maj_allele:[2,1,1,""],min_allele:[2,1,1,""],non_miss:[2,1,1,""],p_mvtest:[2,1,1,""],p_variance:[2,1,1,""],ph_label:[2,1,1,""],pos:[2,1,1,""],print_header:[2,2,1,""],print_result:[2,2,1,""],rsid:[2,1,1,""],runtime:[2,1,1,""],stringify:[2,2,1,""]},"meanvar.mvstandardizer":{Standardizer:[2,3,1,""]},"meanvar.mvstandardizer.Standardizer":{destandardize:[2,2,1,""],standardize:[2,2,1,""]},"meanvar.simple_timer":{SimpleTimer:[2,3,1,""]},"meanvar.simple_timer.SimpleTimer":{report:[2,2,1,""],reset:[2,2,1,""],result:[2,2,1,""],runtime:[2,2,1,""]},"pygwas.bed_parser":{Parser:[3,3,1,""]},"pygwas.bed_parser.Parser":{ReportConfiguration:[3,2,1,""],alleles:[3,1,1,""],bed_file:[3,1,1,""],bim_file:[3,1,1,""],extract_genotypes:[3,2,1,""],fam_file:[3,1,1,""],families:[3,1,1,""],filter_missing:[3,2,1,""],geno_conversions:[3,1,1,""],genotype_file:[3,1,1,""],ind_count:[3,1,1,""],ind_mask:[3,1,1,""],init_genotype_file:[3,2,1,""],load_bim:[3,2,1,""],load_fam:[3,2,1,""],load_genotypes:[3,2,1,""],markers:[3,1,1,""],populate_iteration:[3,2,1,""]},"pygwas.boundary":{BoundaryCheck:[3,3,1,""]},"pygwas.boundary.BoundaryCheck":{LoadExclusions:[3,2,1,""],LoadSNPs:[3,2,1,""],NoExclusions:[3,2,1,""],ReportConfiguration:[3,2,1,""],TestBoundary:[3,2,1,""],beyond_upper_bound:[3,1,1,""],bounds:[3,1,1,""],chrom:[3,1,1,""],dropped_snps:[3,1,1,""],ignored_rs:[3,1,1,""],target_rs:[3,1,1,""],valid:[3,1,1,""]},"pygwas.data_parser":{DataParser:[3,3,1,""],check_inclusions:[3,4,1,""]},"pygwas.data_parser.DataParser":{boundary:[3,1,1,""],compressed_pedigree:[3,1,1,""],get_effa_freq:[3,2,1,""],get_loci:[3,2,1,""],has_fid:[3,1,1,""],has_liability:[3,1,1,""],has_parents:[3,1,1,""],has_pheno:[3,1,1,""],has_sex:[3,1,1,""],ind_exclusions:[3,1,1,""],ind_inclusions:[3,1,1,""],ind_miss_tol:[3,1,1,""],max_maf:[3,1,1,""],min_maf:[3,1,1,""],missing_representation:[3,1,1,""],missing_storage:[3,1,1,""],snp_miss_tol:[3,1,1,""],valid_indid:[3,6,1,""]},"pygwas.exceptions":{InvalidBoundarySpec:[3,5,1,""],InvalidSelection:[3,5,1,""],InvariantVar:[3,5,1,""],MalformedInputFile:[3,5,1,""],NanInResult:[3,5,1,""],NoMatchedPhenoCovars:[3,5,1,""],ReportableException:[3,5,1,""],TooFewAlleles:[3,5,1,""],TooManyAlleles:[3,5,1,""],UnsolvedLocus:[3,5,1,""]},"pygwas.exceptions.TooManyAlleles":{alleles:[3,1,1,""],chr:[3,1,1,""],index:[3,1,1,""],pos:[3,1,1,""],rsid:[3,1,1,""]},"pygwas.impute_parser":{Encoding:[3,3,1,""],Parser:[3,3,1,""],SetEncoding:[3,4,1,""]},"pygwas.impute_parser.Encoding":{Additive:[3,1,1,""],Dominant:[3,1,1,""],Genotype:[3,1,1,""],Raw:[3,1,1,""],Recessive:[3,1,1,""]},"pygwas.impute_parser.Parser":{ReportConfiguration:[3,2,1,""],archives:[3,1,1,""],chroms:[3,1,1,""],current_chrom:[3,1,1,""],current_file:[3,1,1,""],current_info:[3,1,1,""],fam_details:[3,1,1,""],gen_ext:[3,1,1,""],get_effa_freq:[3,2,1,""],get_next_line:[3,2,1,""],info_ext:[3,1,1,""],info_files:[3,1,1,""],info_threshold:[3,1,1,""],load_family_details:[3,2,1,""],load_genotypes:[3,2,1,""],populate_iteration:[3,2,1,""]},"pygwas.locus":{Locus:[3,3,1,""]},"pygwas.locus.Locus":{alleles:[3,1,1,""],chr:[3,1,1,""],exp_hetero_freq:[3,1,1,""],flip:[3,2,1,""],hetero_count:[3,1,1,""],hetero_freq:[3,1,1,""],maf:[3,1,1,""],maj_allele_count:[3,1,1,""],major_allele:[3,1,1,""],min_allele_count:[3,1,1,""],minor_allele:[3,1,1,""],missing_allele_count:[3,1,1,""],p:[3,1,1,""],pos:[3,1,1,""],q:[3,1,1,""],rsid:[3,1,1,""],sample_size:[3,1,1,""],total_allele_count:[3,1,1,""]},"pygwas.mach_parser":{Encoding:[3,3,1,""],Parser:[3,3,1,""]},"pygwas.mach_parser.Encoding":{Dosage:[3,1,1,""]},"pygwas.mach_parser.Parser":{ReportConfiguration:[3,2,1,""],chrpos_encoding:[3,1,1,""],chunk_stride:[3,1,1,""],dosage_ext:[3,1,1,""],get_effa_freq:[3,2,1,""],info_ext:[3,1,1,""],load_family_details:[3,2,1,""],load_genotypes:[3,2,1,""],min_rsquared:[3,1,1,""],openfile:[3,2,1,""],parse_genotypes:[3,2,1,""],populate_iteration:[3,2,1,""]},"pygwas.parsed_locus":{ParsedLocus:[3,3,1,""]},"pygwas.parsed_locus.ParsedLocus":{cur_idx:[3,1,1,""],genotype_data:[3,1,1,""],next:[3,2,1,""]},"pygwas.pedigree_parser":{Parser:[3,3,1,""]},"pygwas.pedigree_parser.Parser":{ReportConfiguration:[3,2,1,""],alleles:[3,1,1,""],datasource:[3,1,1,""],genotypes:[3,1,1,""],get_loci:[3,2,1,""],individual_mask:[3,1,1,""],invalid_loci:[3,1,1,""],load_genotypes:[3,2,1,""],load_mapfile:[3,2,1,""],locus_count:[3,1,1,""],mapfile:[3,1,1,""],markers:[3,1,1,""],markers_maf:[3,1,1,""],populate_iteration:[3,2,1,""],rsids:[3,1,1,""]},"pygwas.pheno_covar":{PhenoCovar:[3,3,1,""]},"pygwas.pheno_covar.PhenoCovar":{add_subject:[3,2,1,""],covariate_data:[3,1,1,""],covariate_labels:[3,1,1,""],destandardize_variables:[3,2,1,""],do_standardize_variables:[3,1,1,""],freeze_subjects:[3,2,1,""],individual_mask:[3,1,1,""],load_covarfile:[3,2,1,""],load_phenofile:[3,2,1,""],missing_encoding:[3,1,1,""],pedigree_data:[3,1,1,""],phenotype_data:[3,1,1,""],phenotype_names:[3,1,1,""],prep_testvars:[3,2,1,""],sex_as_covariate:[3,1,1,""],test_variables:[3,1,1,""]},"pygwas.snp_boundary_check":{SnpBoundaryCheck:[3,3,1,""]},"pygwas.snp_boundary_check.SnpBoundaryCheck":{NoExclusions:[3,2,1,""],ReportConfiguration:[3,2,1,""],TestBoundary:[3,2,1,""]},"pygwas.standardizer":{NoStandardization:[3,3,1,""],StandardizedVariable:[3,3,1,""],get_standardizer:[3,4,1,""],set_standardizer:[3,4,1,""]},"pygwas.standardizer.NoStandardization":{destandardize:[3,2,1,""],standardize:[3,2,1,""]},"pygwas.standardizer.StandardizedVariable":{covar_count:[3,1,1,""],covariates:[3,1,1,""],datasource:[3,1,1,""],destandardize:[3,2,1,""],get_covariate_name:[3,2,1,""],get_covariate_names:[3,2,1,""],get_phenotype_name:[3,2,1,""],get_variables:[3,2,1,""],idx:[3,1,1,""],missing:[3,1,1,""],pheno_count:[3,1,1,""],phenotypes:[3,1,1,""],standardize:[3,2,1,""]},"pygwas.transposed_pedigree_parser":{Parser:[3,3,1,""]},"pygwas.transposed_pedigree_parser.Parser":{ReportConfiguration:[3,2,1,""],filter_missing:[3,2,1,""],load_genotypes:[3,2,1,""],load_tfam:[3,2,1,""],populate_iteration:[3,2,1,""],process_genotypes:[3,2,1,""]},meanvar:{mv_esteq:[2,0,0,"-"],mvresult:[2,0,0,"-"],mvstandardizer:[2,0,0,"-"],simple_timer:[2,0,0,"-"]},pygwas:{BuildReportLine:[3,4,1,""],Exit:[3,4,1,""],ExitIf:[3,4,1,""],bed_parser:[3,0,0,"-"],boundary:[3,0,0,"-"],data_parser:[3,0,0,"-"],exceptions:[3,0,0,"-"],impute_parser:[3,0,0,"-"],locus:[3,0,0,"-"],mach_parser:[3,0,0,"-"],parsed_locus:[3,0,0,"-"],pedigree_parser:[3,0,0,"-"],pheno_covar:[3,0,0,"-"],snp_boundary_check:[3,0,0,"-"],standardizer:[3,0,0,"-"],sys_call:[3,4,1,""],transposed_pedigree_parser:[3,0,0,"-"]}},objnames:{"0":["py","module","Python module"],"1":["py","attribute","Python attribute"],"2":["py","method","Python method"],"3":["py","class","Python class"],"4":["py","function","Python function"],"5":["py","exception","Python exception"],"6":["py","staticmethod","Python static method"]},objtypes:{"0":"py:module","1":"py:attribute","2":"py:method","3":"py:class","4":"py:function","5":"py:exception","6":"py:staticmethod"},terms:{"3rd":3,"abstract":[2,3],"byte":3,"case":[2,3],"class":[2,3],"default":3,"final":3,"function":[1,2,3],"import":2,"int":3,"new":3,"public":[1,3],"return":[2,3],"short":3,"static":3,"switch":3,"throw":3,"true":[2,3],"try":3,"var":[2,3],access:3,accompani:3,accord:3,act:[2,3],actual:[2,3],add:3,add_subject:3,addit:3,after:3,ahead:3,alert:2,algorithm:[2,3],all:[2,3],allel:[2,3],allele2:3,allelel:3,allow:[2,3],along:3,alreadi:3,also:3,amount:3,analysi:[1,2,3],analyz:3,ani:[2,3],anyth:2,appear:3,appli:3,applic:[1,2,3],approach:3,appropri:3,archiv:3,archive_list:3,argument:3,around:3,arrai:3,assign:3,associ:[1,2,3],assum:3,attribut:3,author:[2,3],automat:3,avail:3,avoid:3,back:[2,3],backward:3,base:[2,3],basic:[2,3],becaus:3,bed:3,bed_fil:3,bed_pars:[],been:3,befor:3,begin:3,below:1,ben:3,beta:[2,3],beta_pvalu:2,beta_stderr:2,beta_valu:2,between:3,beyond_upper_bound:3,bim:3,bim_fil:3,binari:3,blin:3,both:3,bound:3,boundari:[],boundarycheck:3,build:3,buildreportlin:3,bunch:3,bvar:3,calcul:[2,3],call:3,can:[2,3],captur:3,chang:3,charact:3,check:3,check_inclus:3,chr:[2,3],chrom:3,chromosom:[2,3],chrpo:3,chrpos_encod:3,chunk:3,chunk_strid:3,circuit:3,client:3,cmd:3,code:3,column:3,combin:2,command:3,comment:3,common:3,complet:[2,3],compon:[2,3],compress:3,compressed_pedigre:3,comput:3,condit:3,configur:3,consid:3,consist:3,consumpt:3,contain:[2,3],content:[],control:3,convent:3,converg:2,convers:3,convert:3,copi:[2,3],correctli:3,could:3,count:[2,3],cov1:[2,3],cov2:3,cov:3,covar:[2,3],covar_count:3,covar_label:2,covari:[2,3],covariate_data:3,covariate_label:3,covn:2,criterion:[2,3],cur_idx:3,current:[2,3],current_chrom:3,current_fil:3,current_info:3,data:[1,2,3],data_pars:[],datapars:3,dataset:[2,3],datasourc:[2,3],defin:3,definit:3,delimit:3,depend:2,desir:3,destandard:[2,3],destandardize_vari:3,destin:3,detail:3,determin:3,didn:3,differ:3,displai:3,dist:3,distanc:3,distinct:3,do_exit:3,do_reset:2,do_standardize_vari:3,doe:3,doesn:3,domin:3,don:3,done:3,dosag:[2,3],dosage_ext:3,dose:3,drop:3,dropped_snp:3,due:3,dure:3,each:[1,2,3],easili:3,eff_alcount:2,effect:[2,3],either:3,ellaps:2,email:3,empti:3,encod:3,encount:3,end:3,end_bound:3,end_r:3,entir:3,entireti:3,enumer:3,equest:2,err:3,errbeta:3,error:[2,3],estim:[2,3],etc:[2,3],evalu:[2,3],even:3,everi:3,everyth:3,exce:3,exceed:3,except:[],exclud:3,exclus:3,execept:3,execut:[2,3],exhaust:3,exit:3,exitif:3,exp_hetero_freq:3,expect:3,explicit:3,extend:1,extens:3,extern:3,extra:3,extract:[1,3],extract_genotyp:3,facilit:[2,3],fail:2,fals:[2,3],fam:3,fam_detail:3,fam_fil:3,famid:3,famili:3,feb:3,few:3,fid:3,file:[1,2,3],filenam:3,filter:3,filter_miss:3,filtrat:3,first:3,fix:3,flag:3,flip:3,follow:[1,2,3],forc:3,form:3,format:3,found:3,fraction:3,fragment:3,freeze_subject:3,frequenc:[2,3],from:[1,2,3],futur:3,gen:3,gen_ext:3,gener:3,genet:3,geno:[2,3],geno_convers:3,genom:3,genotyp:[2,3],genotype_data:3,genotype_fil:3,get_covariate_nam:3,get_effa_freq:3,get_loci:3,get_next_lin:3,get_phenotype_nam:3,get_standard:3,get_vari:3,gnu:1,gpl:1,group:3,gwa:[1,2,3],gzip:3,hand:3,has_fid:3,has_liabl:3,has_par:3,has_pheno:3,has_sex:3,have:3,header:[2,3],here:3,hetero_count:3,hetero_freq:3,heterozygot:3,hold:2,howev:[2,3],html:1,http:1,identifi:3,idx:3,ignor:3,ignored_r:3,iid:3,impli:3,imput:3,impute_pars:[],includ:[1,3],inclus:3,ind_count:3,ind_exclus:3,ind_id:3,ind_inclus:3,ind_mask:3,ind_miss_tol:3,index:3,indic:3,indid:3,individu:3,individual_mask:3,info:3,info_ext:3,info_fil:3,info_threshold:3,inform:3,init_genotype_fil:3,initi:3,input:[2,3],insid:3,instanc:3,instanti:3,intend:3,intens:3,interest:[1,3],interfac:[1,3],intern:3,interpret:3,intial:3,invalid:3,invalid_loci:3,invalidboundaryspec:3,invalidselect:3,invariantvar:3,isn:3,issu:3,item:3,iter:[2,3],itself:3,jan:3,just:3,keep:3,kei:3,keyword:3,kilobase:3,kind:3,know:3,kwarg:[2,3],label:[2,3],larger:3,legaci:3,liabil:3,librari:[],licens:1,like:3,limit:3,line:3,list:[2,3],littl:3,lmpv:2,load:3,load_bim:3,load_covarfil:3,load_fam:3,load_family_detail:3,load_genotyp:3,load_mapfil:3,load_phenofil:3,load_tfam:3,loadexclus:3,loadsnp:3,local:[2,3],locat:3,loci:[2,3],locu:[],locus_count:3,log:[2,3],longer:3,lower:3,mach:3,mach_pars:[],made:3,maf:[2,3],maj:2,maj_allel:2,maj_allele_count:3,major:[2,3],major_allel:3,make:3,malform:3,malformed_boundari:3,malformedinputfil:3,mani:[2,3],map3:3,mapfil:3,marker:3,markers_maf:3,mask:3,massiv:3,master:[],match:3,matrix:3,max_maf:3,mean:[2,3],meaning:3,meaningless:3,meanvar:[],meanvaresteq:2,megabas:3,member:3,memoeri:3,memori:3,messag:[2,3],method:[2,3],mid:3,might:2,min:2,min_allel:2,min_allele_count:3,min_maf:3,min_rsquar:3,minimac:3,minor:[2,3],minor_allel:3,miss:[2,3],missing:3,missing_allele_count:3,missing_encod:3,missing_in_geno:3,missing_represent:3,missing_storag:3,misspel:3,mode:2,model:3,moder:3,modul:[],more:3,most:2,mostli:3,move:3,msg:[2,3],much:3,multipl:3,must:3,mv_esteq:[],mvresult:[],mvstandard:[],name:3,nan:3,naninresult:3,natur:3,necessari:3,next:3,noexclus:3,nomatchedphenocovar:3,non:[2,3],non_miss:2,non_miss_count:2,none:[2,3],nonmiss:3,normal:3,nostandard:3,number:[2,3],numer:3,numpi:3,object:[2,3],observ:3,obtain:3,off:3,offset:3,onc:3,onli:[1,2,3],open:[2,3],openfil:3,optim:2,option:[2,3],order:[2,3],org:1,orient:3,origin:[2,3],other:3,otherwis:3,out:3,outcom:2,output:[2,3],over:3,overrid:3,own:[1,3],p_mvtest:2,p_varianc:2,packag:[],pair:3,parallel:[],param:3,paramet:[2,3],parent:3,pars:3,parse_genotyp:3,parsed_locu:[],parsedlocu:3,parselocu:3,parser:[2,3],part:2,particip:3,pass:3,ped:3,pedigre:3,pedigree_data:3,pedigree_pars:[],per:3,perform:[2,3],permit:3,ph1:3,ph2:3,ph_label:2,pheno:[2,3],pheno_count:3,pheno_covar:[],phenocovar:3,phenofil:3,phenotyp:[2,3],phenotype_data:3,phenotype_nam:3,place:3,placehold:3,plai:3,plink:3,plugin:[2,3],point:3,popul:3,populate_iter:3,portion:3,pos1:3,pos2:3,posit:[2,3],posn:3,possibl:[2,3],pour:3,prefix:3,prep:3,prep_testvar:3,prepar:3,present:3,print:2,print_head:2,print_result:2,prior:3,probabl:3,problem:2,process:3,process_genotyp:3,produc:3,program:1,propag:3,properli:[2,3],provid:[2,3],pull:3,purpos:3,put:3,pvalu:2,pygwa:[],python:3,qualiti:3,ram:3,rang:[2,3],raw:3,reach:3,read:3,readi:3,real:3,realli:3,reason:[],recess:3,record:3,reduc:3,refer:3,regardless:3,region:3,relat:3,releas:1,relev:3,remain:3,remov:3,replac:[2,3],report:[2,3],reportableexcept:3,reportconfigur:3,repres:[1,3],represent:3,requir:[2,3],rescal:[2,3],reset:[2,3],resid:3,result:[2,3],revert:2,right:3,rsid:[2,3],rsquar:3,rule:3,run:[2,3],runanalysi:2,runmeanvar:2,runtim:2,same:3,sampl:3,sample_fil:3,sample_s:3,scale:[2,3],second:[2,3],seen:3,sensic:3,separ:3,set:3,set_standard:3,setencod:3,setup:2,sever:3,sex:3,sex_as_covari:3,should:3,similar:3,simpl:[2,3],simple_tim:[],simpler:3,simpletim:2,simpli:3,sinc:3,singl:[1,2,3],singular:3,situat:[2,3],size:3,snp:[2,3],snp_boundary_check:[],snp_miss_tol:3,snpboundari:3,snpboundarycheck:3,some:[2,3],sort:3,sourc:[2,3],space:3,specfi:3,specif:[2,3],specifi:3,standard:[],standardizedvari:[2,3],start:[2,3],start_bound:3,std:[2,3],stderr:3,stderror:2,stdout:[2,3],still:2,stopiter:3,store:3,stream:2,stride:3,string:[2,3],stringifi:2,strip:2,structur:3,stub:3,studi:3,style:3,subject:3,suffix:3,support:[1,3],sure:3,sval:3,sys_cal:3,tab:3,target:3,target_r:3,task:[],test:3,test_vari:3,testboundari:3,text:3,tfam:3,tfam_fil:3,than:3,them:3,thi:[2,3],thing:3,think:[],those:[1,3],though:3,threshold:3,thrown:3,time:[2,3],todd:3,tol:2,too:3,toofewallel:3,took:2,toomanyallel:3,total:[2,3],total_allele_count:3,total_locus_count:3,tped:3,tped_fil:3,translat:3,transpos:3,transposed_pedigree_pars:[],travers:3,truth:3,turn:3,two:3,type:[1,3],under:1,unlik:3,unsolvedlocu:3,until:3,updat:3,upon:3,upper:3,user:[2,3],util:1,valid:[2,3],valid_indid:3,valu:[2,3],var1:3,var2:3,variabl:[2,3],varianc:2,variant:3,variou:3,verbos:2,version:[1,3],via:3,wai:[2,3],well:3,were:[2,3],what:3,when:[2,3],where:3,which:[1,3],whitespac:3,who:[1,3],wish:1,within:[2,3],without:3,won:3,work:[2,3],would:3,write:[2,3],www:1,you:3,zhang:3},titles:["MVtest API","MVtest API Libraries","meanvar package","pygwas package"],titleterms:{api:[0,1],bed_pars:3,boundari:3,content:[2,3],data_pars:3,except:3,impute_pars:3,librari:1,locu:3,mach_pars:3,meanvar:2,modul:[2,3],mv_esteq:2,mvresult:2,mvstandard:2,mvtest:[0,1],packag:[2,3],parsed_locu:3,pedigree_pars:3,pheno_covar:3,pygwa:3,simple_tim:2,snp_boundary_check:3,standard:3,transposed_pedigree_pars:3}})