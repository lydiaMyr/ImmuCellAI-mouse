
#' Title
#'
#' @param sample
#' @param data_type
#' @param group_tag
#' @param customer
#' @param sig_file
#' @param exp_file
#'
#' @return
#' @export
#'
#' @examples
ImmuCellAI_mouse<-function(sample,data_type,group_tag,customer,sig_file=NULL,exp_file=NULL){
  data("l1_marker")
  data("l2_marker")
  data("l3_marker")
  data("l1_cell_correction_matrix_new")
  data("l2_cell_correction_matrix_new")
  data("l3_cell_correction_matrix_new")
  group_index=0
  if (group_tag){
    #group_index=as.numeric(as.vector(unlist(grep("group",row.names(sample)))))
    group_column<-sample[1,]
    group_content<<-sample[1,]
    sample=sample[-1,]
  }
  sam = apply(sample,2,as.numeric)
  row.names(sam) = row.names(sample)
  # tt = intersect(row.names(sam),as.vector(unlist(paper_marker)))
  # genes = intersect(tt,row.names(marker_exp))
  sam_exp = as.matrix(sam)
  colnames(sam_exp) = colnames(sample)
  sample_exp<-sam_exp
  paper_marker<-l1_marker

  ref_pre<-gsva(as.matrix(marker_exp_raw[,names(paper_marker)]),paper_marker,method="ssgsea",ssgsea.norm=TRUE)
  cell_cor_new_mat<-matrix(rep(0,length(names(l1_marker))*length(names(l1_marker))),ncol=length(names(l1_marker)))
  row.names(cell_cor_new_mat)=names(l1_marker)
  colnames(cell_cor_new_mat)=names(l1_marker)
  diag(cell_cor_new_mat)<-1
  sample=as.matrix(sample_exp)
  sam = apply(sample,2,as.numeric)
  row.names(sam) = row.names(sample_exp)
  tt = intersect(row.names(sam),as.vector(unlist(paper_marker)))
  genes = intersect(tt,row.names(marker_exp))
  sam_exp = as.matrix(sam[genes,])
  row.names(sam_exp) = genes
  colnames(sam_exp) = colnames(sample)
  # marker_exp = marker_exp[genes,]
  tt = intersect(row.names(sam_exp),as.vector(unlist(paper_marker)))
  genes = intersect(tt,row.names(marker_exp))
  marker_exp_partial<-marker_exp[genes,]
  marker_tag_mat = c()
  # binary tag matrix
  for(cell in names(paper_marker)){
    tag = marker_tag(genes,as.vector(unlist(paper_marker[cell])))
    marker_tag_mat = cbind(marker_tag_mat,tag)
  }
  row.names(marker_tag_mat) = row.names(marker_exp_partial)
  colnames(marker_tag_mat) = names(paper_marker)
  marker_reccur<-marker_times_count(paper_marker)
  gene_weight_mt<-gene_ref_weight_mt(paper_marker,marker_reccur)
  gene_weight_mt<-gene_weight_mt[genes,]
  immune_deviation_sample = apply(sam_exp,2,function(x) sample_ratio(x,marker_exp_partial[,names(paper_marker)],marker_tag_mat,data_type,gene_weight_mt,ref_pre,cell_cor_new_mat))
  #row.names(immune_deviation_sample)<-row.names(marker_exp)
  infil_marker<-genes
  # infil_exp<-data.frame(as.vector(unlist(immune_deviation_sample[genes,])))
  # row.names(infil_exp)<-as.vector(unlist(sapply(colnames(immune_deviation_sample),function(x) paste(x,genes,sep="_"))))
  # colnames(infil_exp)<-"Ratio"
  # infil_marker<-list()
  # for(sam in colnames(immune_deviation_sample)){
  #   infil_marker[[sam]]<-row.names(infil_exp)[grep(sam,row.names(infil_exp))]
  # }
  # infil_exp$Ratio<-as.numeric(infil_exp$Ratio)
  infil_mt<-c()
  for(cell in names(l1_marker)[-c(5,6)]){
    genes_new<-intersect(l1_marker[[cell]],row.names(immune_deviation_sample))
    infil_marker_cell<-list()
    infil_exp<-data.frame(as.vector(unlist(immune_deviation_sample[genes_new,])))
    row.names(infil_exp)<-as.vector(unlist(sapply(colnames(immune_deviation_sample),function(x) paste(x,genes_new,sep="_"))))
    colnames(infil_exp)<-"Ratio"
    for(sam in colnames(immune_deviation_sample)){
      infil_marker_cell[[sam]]<-row.names(infil_exp)[grep(sam,row.names(infil_exp))]
    }
    infil_score<-gsva(as.matrix(infil_exp),infil_marker_cell,method="ssgsea",ssgsea.norm=TRUE)
    infil_mt<-cbind(infil_mt,infil_score)
  }
  infil_mt[which(infil_mt<0)]=0
  infil_score<-apply(infil_mt,1,function(x) sum(x))
  #infil_score[which(infil_score>2)]=runif(length(which(infil_score>2)),1,2)
  #k=(0.9-0.1)/(max(infil_score)-min(infil_score))
  #infil_score<-0.1+k*(infil_score-min(infil_score))


  tt = intersect(row.names(sample_exp),as.vector(unlist(paper_marker)))
  genes = intersect(tt,row.names(marker_exp))
  sam_exp<-sample_exp[genes,]
  # layer1_result<-getResult(sam_exp,data_type,marker_exp_raw[,names(paper_marker)],paper_marker,ref_pre,cell_cor_new_mat)
  # layer1_result<-round(t(apply(layer1_result,1,function(x) x/sum(x))),4)
  layer1_pre<-getResult(sam_exp,data_type,marker_exp_raw[,names(paper_marker)],paper_marker,ref_pre,cell_cor_new_mat,l1_cell_correction_matrix_new)
  #infiltration_score<-apply(layer1_result,1,sum)
  layer1_result<-round(t(apply(layer1_pre,1,function(x) x/sum(x))),4)

  marker_exp<-marker_exp_raw
  paper_marker<-l2_marker
  cell_cor_mat = matrix(rep(0,length(names(paper_marker))*length(names(paper_marker))),ncol=length(names(paper_marker)))
  row.names(cell_cor_mat) = names(paper_marker)
  colnames(cell_cor_mat) = names(paper_marker)
  # binary tag matrix
  #print(colnames(marker_exp_raw))
  #print(names(paper_marker))
  ref_pre<-gsva(as.matrix(marker_exp_raw[,names(paper_marker)]),paper_marker,method="ssgsea",ssgsea.norm=TRUE)
  tt = intersect(row.names(sample_exp),as.vector(unlist(paper_marker)))
  genes = intersect(tt,row.names(marker_exp))
  sam_exp<-sample_exp[genes,]
  # layer2_pre<-getResult(sam_exp,data_type,marker_exp_raw[,names(paper_marker)],paper_marker,ref_pre,cell_cor_new_mat)
  # layer2_pre_trans<-data.frame(t(layer2_pre),check.names=F)
  layer2_pre<-getResult(sam_exp,data_type,marker_exp_raw[,names(paper_marker)],paper_marker,ref_pre,cell_cor_new_mat,l2_cell_correction_matrix_new)
  # layer2_pre1<-t(getResult1(sample_exp,data_type,marker_exp_raw[,names(paper_marker)],paper_marker,ref_pre,cell_cor_new_mat,"l2"))

  layer2_pre_trans<-data.frame(t(layer2_pre),check.names=F)


  tt<-t(layer_norm(layer2_pre_trans,T_sub_ls))
  tmp<-round(t(apply(tt,1,function(x) x/sum(x))),4)
  layer2_T_norm<-tmp*layer1_result[,"T_cell"]
  tt<-t(layer_norm(layer2_pre_trans,B_sub_ls))
  tmp<-round(t(apply(tt,1,function(x) x/sum(x))),4)
  layer2_B_norm<-tmp*layer1_result[,"B_cell"]
  tt<-t(layer_norm(layer2_pre_trans,DC_sub_ls))
  tmp<-round(t(apply(tt,1,function(x) x/sum(x))),4)
  layer2_DC_norm<-tmp*layer1_result[,"Dendritic_cells"]
  tt<-t(layer_norm(layer2_pre_trans,Gra_sub_ls))
  tmp<-round(t(apply(tt,1,function(x) x/sum(x))),4)
  layer2_Gra_norm<-tmp*layer1_result[,"Granulocytes"]
  tt<-t(layer_norm(layer2_pre_trans,Macro_sub_ls))
  tmp<-round(t(apply(tt,1,function(x) x/sum(x))),4)
  layer2_macro_norm<-tmp*layer1_result[,"Macrophage"]
  #
  # layer2_T_norm<-t(layer_norm(layer2_pre_trans,names(T_sub_ls)))*layer1_result[,"T_cell"]
  # layer2_B_norm<-t(layer_norm(layer2_pre_trans,names(B_sub_ls)))*layer1_result[,"B_cell"]
  # layer2_DC_norm<-t(layer_norm(layer2_pre_trans,names(DC_sub_ls)))*layer1_result[,"Dendritic_cells"]
  # layer2_Gra_norm<-t(layer_norm(layer2_pre_trans,names(Gra_sub_ls)))*layer1_result[,"Granulocytes"]
  # layer2_macro_norm<-t(layer_norm(layer2_pre_trans,names(Macro_sub_ls)))*layer1_result[,"Macrophage"]
  #layer3
  #Tcm T helper鐨刴arker闇€瑕佽皟鏁?
  layer2_cell_all<-cbind(layer2_T_norm,layer2_B_norm,layer2_DC_norm,layer2_Gra_norm,layer2_macro_norm)

  marker_exp<-marker_exp_raw
  paper_marker<-l3_marker



  ref_pre<-gsva(as.matrix(marker_exp_raw[,names(paper_marker)]),paper_marker,method="ssgsea",ssgsea.norm=TRUE)
  tt = intersect(row.names(sample_exp),as.vector(unlist(paper_marker)))
  genes = intersect(tt,row.names(marker_exp))
  sam_exp<-sample_exp[genes,]
  # layer3_pre<-getResult(sam_exp,data_type,marker_exp_raw[,names(paper_marker)],paper_marker,ref_pre,cell_cor_new_mat)
  # layer3_pre_trans<-data.frame(t(layer3_pre),check.names=F)
  layer3_pre<-getResult(sam_exp,data_type,marker_exp_raw[,names(paper_marker)],paper_marker,ref_pre,cell_cor_new_mat,l3_cell_correction_matrix_new)
  #
  layer3_pre_trans<-data.frame(t(layer3_pre),check.names=F)

  # layer3_CD4_sub<-t(layer_norm(layer3_pre_trans,names(CD4_T_sub_ls)))
  # # tmp<-round(t(apply(tt,1,function(x) x/sum(x))),4)
  # # layer3_CD4_sub<-tmp*layer2_T_norm[,"CD4_T_cell"]
  #
  # layer3_CD8_sub<-t(layer_norm(layer3_pre_trans,names(CD8_T_sub_ls)))
  # # tmp<-round(t(apply(tt,1,function(x) x/sum(x))),4)
  # # layer3_CD8_sub<-tmp*layer2_T_norm[,"CD8_T_cell"]

  tt<-t(layer_norm(layer3_pre_trans,CD4_T_sub_ls))
  tmp<-round(t(apply(tt,1,function(x) x/sum(x))),4)
  layer3_CD4_sub<-tmp*layer2_T_norm[,"CD4_T_cell"]

  tt<-t(layer_norm(layer3_pre_trans,CD8_T_sub_ls))
  tmp<-round(t(apply(tt,1,function(x) x/sum(x))),4)
  layer3_CD8_sub<-tmp*layer2_T_norm[,"CD8_T_cell"]

  #result_final<-cbind(layer1_result,layer2_T_norm,layer2_B_norm,layer2_DC_norm,layer2_Gra_norm,layer2_macro_norm,layer3_CD4_sub,layer3_CD8_sub)
  result_final<-cbind(layer1_result,layer2_cell_all,layer3_CD4_sub,layer3_CD8_sub)
  result_final<-round(result_final,4)
  result_final[grep(TRUE,is.na(result_final))]=0
  result_mat<-result_final
  if(group_tag){
    group_name<-sort(unique(as.vector(unlist(group_content))))
    p_va=c()
    group_column<-as.numeric(as.factor(as.vector(unlist(group_content))))
    Infiltration_score<-round(infil_score,3)
    result_group=cbind(result_mat,Infiltration_score,group_column)
    result_tt=apply(result_group,2,as.numeric)
    # print("done")
    if (length(group_name)>2){
      for (cell in colnames(result_group)){
        result_group_new<-result_group[,c(cell,"group_column")]
        t=aov(group_column~.,data.frame(result_group_new))
        p_va=c(p_va,round(summary(t)[[1]][["Pr(>F)"]],2))
      }
    }else{
      g1_index=grep(1,group_column)
      g2_index=grep(2,group_column)
      result_mat1<-cbind(result_mat,Infiltration_score)
      for (cell in colnames(result_mat)){
        c_=wilcox.test(result_mat[g1_index,cell],result_mat[g2_index,cell])
        p_va=c(p_va,round(c_$p.value,2))
      }
    }
    #  print("done1")
    #p_va=p_va[-26]
    row.names(result_tt)=row.names(result_group)
    result_tt=data.frame(result_tt)
    #print(res)
    exp_median=aggregate(.~group_column,data=result_tt,median)

    exp_median=rbind(exp_median[,-1],p_va)
    #row.names(exp_median)=c(group_name,"p value")
    row.names(exp_median)=c(group_name,"p value")
    group_fre<<-exp_median
    write.table(group_fre,save_group,sep="\t",quote=F,col.names = NA)
    Infiltration_score<-round(infil_score,3)
    T_FRE<<-cbind(result_mat,Infiltration_score)
    #print("done2")
    plot_fun("boxplot",1,save_plot)
  }
  Infiltration_score<-round(infil_score,3)
  T_FRE<<-cbind(result_mat,Infiltration_score)
  write.table(T_FRE,save_abun,sep="\t",quote=F,col.names = NA)
}

#' Title
#'
#' @param sample
#' @param data_type
#' @param marker_exp
#' @param paper_marker
#' @param ref_pre
#' @param cell_cor_new_mat
#' @param compensation_matrix
#'
#' @return
#' @export
#'
#' @examples
getResult = function(sample,data_type,marker_exp,paper_marker,ref_pre,cell_cor_new_mat,compensation_matrix){
  #marker_exp <- read.csv("/project/xiamx/immunecell/data/GEO/main/prepare_for_train_set.txt",sep = "\t",row.names = 1)
  #marker_exp <-
  sample=as.matrix(sample)
  sam = apply(sample,2,as.numeric)
  row.names(sam) = row.names(sample)
  tt = intersect(row.names(sam),as.vector(unlist(paper_marker)))
  genes = intersect(tt,row.names(marker_exp))
  sam=data.frame(sam)
  sam_exp = as.matrix(sam[genes,])
  colnames(sam_exp) = colnames(sample)
  marker_exp = marker_exp[genes,]


  marker_tag_mat = c()
  # binary tag matrix
  for(cell in names(paper_marker)){
    tag = marker_tag(genes,as.vector(unlist(paper_marker[cell])))
    marker_tag_mat = cbind(marker_tag_mat,tag)
  }
  row.names(marker_tag_mat) = row.names(marker_exp)
  colnames(marker_tag_mat) = names(paper_marker)
  marker_reccur<-marker_times_count(paper_marker)
  # gene_weight_mt<-gene_ref_weight_mt(paper_marker,marker_reccur)
  # gene_weight_mt<-gene_weight_all_ls[[num_layer]]
  # genes<-intersect(genes,row.names(gene_weight_mt))
  gene_weight_mt<-c()

  exp_new = apply(sam_exp,2,function(x) sample_ratio(x,marker_exp[,names(paper_marker)],marker_tag_mat,data_type,gene_weight_mt,ref_pre,cell_cor_new_mat))
  #print(head(exp_new))
  #print(paper_marker)
  result = gsva(exp_new,paper_marker,method="ssgsea",ssgsea.norm=TRUE)
  #exp_raio_new<-apply(sam_exp,2,function(x) sample_ratio_correction(x,marker_exp[,names(paper_marker)],marker_tag_mat,data_type,gene_weight_mt))
  #print(head(result))
  #
  result<-apply(result,1,function(x) round(x-min(x)+min(abs(x)),3))

  # if(ncol(result)<3){
  #   result[which(result<0)]=0
  # }else{
  #   result = result - apply(result,1,min)
  #   #result[which(result<0)]=0
  # }
  compensation_matrix_num = apply(compensation_matrix,2,as.numeric)
  # progress$set(value = 20,detail = "Adjusting result by Compensation matrix")
  # incProgress(0.2, detail = "Immune infiltration calculating")
  # Sys.sleep(0.5)
  # shinyWidgets::updateProgressBar(title = 'Immune infiltration calculating',id = "pb2",value=70,session=getDefaultReactiveDomain())
  # if(customer==0){
  row.names(compensation_matrix_num) = row.names(compensation_matrix)
  result_mat = compensation(t(result),compensation_matrix_num)
  if(ncol(result_mat)>1){
    result_mat=apply(result_mat,1,function(x) round(x,3))
  }else{
    result_mat=t(round(result_mat,3))
  }
  return(result_mat)
}



#' Title
#'
#' @param raw_score
#' @param compensation_matrix
#'
#' @return
#' @export
#'
#' @examples
compensation = function(raw_score,compensation_matrix){
  raw_score=as.matrix(raw_score)
  compensation_matrix = compensation_matrix
  diag(compensation_matrix) = 1
  rows <- rownames(raw_score)[rownames(raw_score) %in%  rownames(compensation_matrix)]
  #print(rows)
  if(ncol(raw_score)==1){
    scores <- as.matrix(pracma::lsqlincon(compensation_matrix[rows,rows], raw_score, lb = 0))
  }else{
    scores <- apply(raw_score[rows,], 2, function(x) pracma::lsqlincon(compensation_matrix[rows,rows], x, lb = 0))

  }
  scores<-apply(scores,1,function(x) round(x-min(x)+min(abs(x)),3))
  #scores[scores < 0] = 0
  colnames(scores) <- rows
  return(t(scores))
}


#' Title
#'
#' @param pre_result
#' @param cell_ls
#'
#' @return
#' @export
#'
#' @examples
layer_norm<-function(pre_result,cell_ls){
  pre_result%>%
    dplyr::mutate(cellType=row.names(.))%>%
    dplyr::filter(cellType%in%cell_ls)->sub_fra
  row.names(sub_fra)<-sub_fra$cellType
  sub_fra<-sub_fra[,-ncol(sub_fra)]
  return(sub_fra)
}


#' Title
#'
#' @param paper_marker
#'
#' @return
#' @export
#'
#' @examples
marker_times_count=function(paper_marker){
  print(test)
  #all_genes=unique(as.vector(unlist(paper_marker)))
  gene_count=list()
  for(cell in names(paper_marker)){
    for(gene in paper_marker[[cell]]){
      gene_count[[gene]]<-c(gene_count[[gene]],cell)
      # if(length(which(names(gene_count)==gene))==0){
      #   gene_count[[gene]]=1
      # }else{
      #   gene_count[[gene]]=gene_count[[gene]]+1
      # }
    }
  }
  return(gene_count)
}

#' Title
#'
#' @param marker_ls
#' @param marker_reccur
#'
#' @return
#' @export
#'
#' @examples
gene_ref_weight_mt<-function(marker_ls,marker_reccur){
  gene_weight_mt<-c()
  for(gene in names(marker_reccur)){
    tmp<-rep(1,length(names(marker_ls)))
    names(tmp)<-names(marker_ls)
    cells<-marker_reccur[[gene]]
    if(length(cells)>1){
      cell_exp_sum<-sum(as.numeric(as.vector(unlist(marker_exp[gene,cells]))))
      for(cell in cells){
        tmp[cell]=marker_exp[gene,cell]/cell_exp_sum
      }
    }
    gene_weight_mt<-rbind(gene_weight_mt,tmp)
  }
  row.names(gene_weight_mt)<-names(marker_reccur)
  colnames(gene_weight_mt)<-names(marker_ls)
  return(gene_weight_mt)
}


#' Title
#'
#' @param marker_exp
#' @param marker_ls
#'
#' @return
#' @export
#'
#' @examples
cell_type_correlation<-function(marker_exp,marker_ls){
  exp_cor_mt<-c()
  for(c1 in names(marker_ls)){
    cor_va<-c()
    c1_exp<-marker_exp[marker_ls[[c1]],c1]
    for(c2 in names(marker_ls)){
      c2_exp<-marker_exp[marker_ls[[c1]],c2]
      cor_tmp<-cor.test(as.vector(unlist(c1_exp)),as.vector(unlist(c2_exp)))
      cor_va<-c(cor_va,cor_tmp$estimate)
    }
    exp_cor_mt<-rbind(exp_cor_mt,cor_va)
  }
  colnames(exp_cor_mt)<-names(marker_ls)
  row.names(exp_cor_mt)<-names(marker_ls)
  return(exp_cor_mt)
}




#' Title
#'
#' @param exp
#' @param gene_name
#' @param data_type
#' @param pre_result
#'
#' @return
#' @export
#'
#' @examples
immune_infiltate_calculate=function(exp,gene_name,data_type,pre_result){
  inf = 0
  names(exp) = gene_name
  for (cell in names(immune_infiltate_marker)){
    abun = 0
    markers = as.vector(unlist(immune_infiltate_marker[cell]))
    for (gene in markers){
      if(data_type == "microarray"){
        abun = abun + as.numeric(exp[gene])/marker_exp_T[gene,cell]
      }else{
        abun = abun + as.numeric(log2(exp[gene]+1))/marker_exp_T[gene,cell]
      }
    }
    inf = inf+abun/length(as.vector(unlist(immune_infiltate_marker[cell])))
  }
  return(inf)
}
#' Title
#'
#' @param comgenes
#' @param tag_gene
#'
#' @return
#' @export
#'
#' @examples
marker_tag = function(comgenes,tag_gene){
  a = comgenes
  a[which(comgenes%in%tag_gene)] = 1
  a[which(a!=1)] = 0
  a = as.numeric(a)
  return(a)
}

#' Title
#'
#' @param data
#' @param marker_exp
#' @param marker_tag_mat
#' @param data_type
#' @param gene_weight_mt
#' @param ref_pre
#' @param cell_cor_new_mat
#'
#' @return
#' @export
#'
#' @examples
sample_ratio = function(data,marker_exp,marker_tag_mat,data_type,gene_weight_mt,ref_pre,cell_cor_new_mat){
  exp = 0
  #print(data[1:2])
  #Eprint(data_type)
  if(data_type == "microarray"){
    for (cell in colnames(marker_exp)){
      # cor_cell<-intersect(names(which(cell_cor_new_mat[cell,]==1)),names(which(ref_pre[,cell]>0)))
      # cor_cell_ratio<-ref_pre[cor_cell,cell]/sum(ref_pre[cor_cell,cell])
      # names(cor_cell_ratio)<-cor_cell
      # for(cc in cor_cell){
      exp<-exp+data/(marker_exp[,cell])*marker_tag_mat[,cell]
    }
    #}
  }else{
    for (cell in colnames(marker_exp)){
      #exp = exp+log2(data+1)*gene_weight_mt[,cell]/marker_exp[,cell]*marker_tag_mat[,cell]
      # cor_cell<-intersect(names(which(cell_cor_new_mat[cell,]==1)),names(which(ref_pre[,cell]>0)))
      # cor_cell_ratio<-ref_pre[cor_cell,cell]/sum(ref_pre[cor_cell,cell])
      # names(cor_cell_ratio)<-cor_cell
      # for(cc in cor_cell){
      exp<-exp+log2(data+1)/(marker_exp[,cell])*marker_tag_mat[,cell]
      #}
    }
  }
  #  print(exp[1:5])
  return(exp)
}
