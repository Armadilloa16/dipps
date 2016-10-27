# Copyright (C) 2016 Lyron Winderbaum
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

# Plots a spatial image
spatialPlot <- function(peaklist_in,fExists_in,
                        plot_var = "intensity",
                        plot_var_transform = "none",
                        plot_var_type = "continuous",
                        mult_peaks = "average",
                        save_plot = FALSE,
                        plot_name_in = "",
                        minX_in = 1,
                        minY_in = 1,
                        display_pixel_borders = FALSE,
                        display_legend = TRUE,
                        return_mI.m = FALSE
){

  if (plot_var_type == "continuous"){
    if (is.na(match(plot_var,c("m.z","SN","QualityFactor","Resolution","intensity","area","count")))){
      print(paste("Plotting Variable",plot_var,"is not currently supported."))
      print("Defaulting to intensity")
      plot_var <- "intensity"
      print("In order to plot non-standard variables the spatialPlot() function will need to be modified.")
    } else {
      if (plot_var == "count" & mult_peaks != "sum") {
        mult_peaks <- "sum"
      }
    }
    if (is.na(match(mult_peaks,c("average","sum","max")))){
      print(paste("Method for reducing multiple peaks",mult_peaks,"is not currently supported."))
      print("Defaulting to averaging")
      mult_peaks <- "average"
      print("In order to use non-standard methods the spatialPlot() function will need to be modified.")
    }
  } else if (plot_var_type != "categorical") {
    print(paste("Variable type",plot_var_type,"is not currently supported."))
    print("Defaulting to categorical")
    plot_var_type <- "categorical"
    print("In order to use non-standard methods the spatialPlot() function will need to be modified.")
  }

  # Deal with multiple peaks.
  temp = as.vector(table(peaklist_in$Acquisition))
  if (plot_var_type == "continuous"){
    if (plot_var == "count") {
      peaklist_in$count = 1
    }
    if (sum(temp>1) > 0) {
      peaklist_in <- switch(mult_peaks,
                            average = switch(plot_var,
                                             m.z           = ddply(peaklist_in,"Acquisition",summarise,m.z = mean(m.z)),
                                             SN            = ddply(peaklist_in,"Acquisition",summarise,SN = mean(SN)),
                                             QualityFactor = ddply(peaklist_in,"Acquisition",summarise,QualityFactor = mean(QualityFactor)),
                                             Resolution    = ddply(peaklist_in,"Acquisition",summarise,Resolution = mean(Resolution)),
                                             intensity     = ddply(peaklist_in,"Acquisition",summarise,intensity = mean(intensity)),
                                             area          = ddply(peaklist_in,"Acquisition",summarise,area = mean(area))
                            ),
                            sum = switch(plot_var,
                                         m.z           = ddply(peaklist_in,"Acquisition",summarise,m.z = sum(m.z)),
                                         SN            = ddply(peaklist_in,"Acquisition",summarise,SN = sum(SN)),
                                         QualityFactor = ddply(peaklist_in,"Acquisition",summarise,QualityFactor = sum(QualityFactor)),
                                         Resolution    = ddply(peaklist_in,"Acquisition",summarise,Resolution = sum(Resolution)),
                                         intensity     = ddply(peaklist_in,"Acquisition",summarise,intensity = sum(intensity)),
                                         area          = ddply(peaklist_in,"Acquisition",summarise,area = sum(area)),
                                         count         = ddply(peaklist_in,"Acquisition",summarise,count = sum(count))
                            ),
                            max = switch(plot_var,
                                         m.z           = ddply(peaklist_in,"Acquisition",summarise,m.z = max(m.z)),
                                         SN            = ddply(peaklist_in,"Acquisition",summarise,SN = max(SN)),
                                         QualityFactor = ddply(peaklist_in,"Acquisition",summarise,QualityFactor = max(QualityFactor)),
                                         Resolution    = ddply(peaklist_in,"Acquisition",summarise,Resolution = max(Resolution)),
                                         intensity     = ddply(peaklist_in,"Acquisition",summarise,intensity = max(intensity)),
                                         area          = ddply(peaklist_in,"Acquisition",summarise,area = max(area))
                            )
      )
    }
  } else if (sum(temp>1) > 0) {
    stop("In order to plot categorical variables spectra must be uniquely specified.")
  }

  if (plot_var_type == "continuous"){
    if (is.na(match(plot_var_transform,c("none","log")))){
      print(paste("Transformation",plot_var_transform,"is not currently supported."))
      print("Defaulting to none")
      plot_var_transform = "none"
      print("In order to use non-standard transformations the spatialPlot() function will need to be modified.")
    }
    if (plot_var_transform == "log"){
      min_plot_var = min(round(log(1 + peaklist_in[,plot_var]),5))
      max_plot_var = max(round(log(1 + peaklist_in[,plot_var]),5))
    } else {
      min_plot_var = min(peaklist_in[,plot_var])
      max_plot_var = max(peaklist_in[,plot_var])
    }
  } else {
    plot_var_transform <- "none"
    if (is.factor(peaklist_in[,plot_var])) {
      peaklist_in[,plot_var] = levels(peaklist_in[,plot_var])[as.numeric(peaklist_in[,plot_var])]
    }
  }

  if (is.data.frame(fExists_in)){
    fExists_in = list(fExists_in)
  }
  mI.m_out <- as.list(1:length(fExists_in))
  for (region_idx in 1:length(fExists_in)){
    fExists = fExists_in[[region_idx]]
    if (length(fExists_in) != 1){
      if (length(fExists_in) == length(plot_name_in)){
        plot_name = plot_name_in[region_idx]
      } else {
        plot_name = toString(region_idx)
      }
      if (length(fExists_in) == length(minX_in)){
        minX = minX_in[region_idx]
      } else {
        minX = 1
      }
      if (length(fExists_in) == length(minY_in)){
        minY = minY_in[region_idx]
      } else {
        minY = 1
      }
    } else {
      plot_name = plot_name_in
      minX = minX_in
      minY = minY_in
    }

    # Generate an image matrix mI
    mI <- fExists

    mI[mI > 0] <- match(mI[mI > 0],peaklist_in$Acquisition)
    subset_of_peaklist = mI[mI > 0 & !is.na(mI)]
    mI[mI > 0 & !is.na(mI)] <- peaklist_in[subset_of_peaklist,plot_var]

    mI$X = (1:nrow(mI)) + minX - 2
    mI.m = melt(mI,id.var="X")
    mI.m$Y = (as.numeric(substring(mI.m$variable,2))) + minY - 2

    if (plot_var_transform == "log"){
      mI.m$value = round(log(1 + mI.m$value),5)
    }

    # Make a new variable empty for acquisition regions that have no peaks,
    #    and for non-acquisition regions set value to NA
    mI.m$empty = FALSE
    if (sum(is.na(mI.m$value)) > 0) {
      mI.m[is.na(mI.m$value),]$empty = TRUE
    }
    if (sum(mI.m[!is.na(mI.m$value),"value"]==0) > 0){
      mI.m[replace(mI.m$value,is.na(mI.m$value),1)==0,]$value = NA
    }

    if (plot_var_type == "categorical"){
      mI.m$value <- factor(mI.m$value)
    }

    if (!return_mI.m){
      p = ggplot(mI.m,aes(X,Y))
      if (display_pixel_borders){
        p = p + geom_tile(data = mI.m,aes(fill=value,alpha=as.numeric(!is.na(value))),colour="grey")
      } else {
        p = p + geom_tile(data = mI.m,aes(fill=value,alpha=as.numeric(!is.na(value))),colour=NA)
      }
      if (display_legend){
        p = p + guides(alpha = FALSE)
      } else {
        p = p + guides(alpha = FALSE,fill=FALSE)
      }
      p = p + geom_tile(data = mI.m,alpha=0.5*as.numeric(mI.m$empty))
      if (plot_var_type == "continuous"){
        p = p + scale_fill_gradient(name=plot_var, limits = c(min_plot_var,max_plot_var))
      }
      #     p = p + ggtitle(plot_name)
      p = p + coord_fixed()
      p = p + scale_y_reverse()

      if (plot_name != ""){
        plot_name = paste(plot_name,"_",sep="")
      }

      if (save_plot){
        ggsave(paste(plot_name,paste(mult_peaks,plot_var,"transformation",plot_var_transform,sep="_"),".png",sep=""),p)
      }

      if (region_idx == 1){
        p_list = list(p)
      } else {
        p_list[[region_idx]] = p
      }
    } else {
      mI.m_out[[region_idx]] = mI.m
    }
  }
  if (return_mI.m){
    if (length(fExists_in)==1){
      return(mI.m_out[[1]])
    } else {
      return(mI.m_out)
    }
  } else {
    if (length(fExists_in)==1){
      return(p_list[[1]])
    } else {
      return(p_list)
    }
  }
}

# Plots an acquisition plot (with acquisition order on the x-axis.)
acquisitionPlot <- function(peaklist_in,
                            plot_var = "intensity",
                            plot_var_transform = "none",
                            mult_peaks = "average",
                            save_plot = FALSE,
                            plot_name = ""){

  if (is.na(match(plot_var,c("m.z","SN","QualityFactor","Resolution","intensity","area","count")))){
    print(paste("Plotting Variable",plot_var,"is not currently supported."))
    print("Defaulting to intensity")
    plot_var <- "intensity"
    print("In order to plot non-standard variables the spatialPlot() function will need to be modified.")
  } else {
    if (plot_var == "count" & mult_peaks != "sum") {
      mult_peaks <- "sum"
    }
  }
  if (is.na(match(mult_peaks,c("average","sum","max")))){
    print(paste("Method for reducing multiple peaks",mult_peaks,"is not currently supported."))
    print("Defaulting to averaging")
    mult_peaks <- "average"
    print("In order to use non-standard methods the spatialPlot() function will need to be modified.")
  }

  # Generate an image matrix mI
  temp = as.vector(table(peaklist_in$Acquisition))
  # Deal with multiple peaks.
  if (plot_var == "count") {
    peaklist_in$count = 1
  }
  if (sum(temp>1) > 0) {
    peaklist_in <- switch(mult_peaks,
                          average = switch(plot_var,
                                           m.z           = ddply(peaklist_in,c("Acquisition","Peaklist"),summarise,m.z = mean(m.z)),
                                           SN            = ddply(peaklist_in,c("Acquisition","Peaklist"),summarise,SN = mean(SN)),
                                           QualityFactor = ddply(peaklist_in,c("Acquisition","Peaklist"),summarise,QualityFactor = mean(QualityFactor)),
                                           Resolution    = ddply(peaklist_in,c("Acquisition","Peaklist"),summarise,Resolution = mean(Resolution)),
                                           intensity     = ddply(peaklist_in,c("Acquisition","Peaklist"),summarise,intensity = mean(intensity)),
                                           area          = ddply(peaklist_in,c("Acquisition","Peaklist"),summarise,area = mean(area))
                          ),
                          sum = switch(plot_var,
                                       m.z           = ddply(peaklist_in,c("Acquisition","Peaklist"),summarise,m.z = sum(m.z)),
                                       SN            = ddply(peaklist_in,c("Acquisition","Peaklist"),summarise,SN = sum(SN)),
                                       QualityFactor = ddply(peaklist_in,c("Acquisition","Peaklist"),summarise,QualityFactor = sum(QualityFactor)),
                                       Resolution    = ddply(peaklist_in,c("Acquisition","Peaklist"),summarise,Resolution = sum(Resolution)),
                                       intensity     = ddply(peaklist_in,c("Acquisition","Peaklist"),summarise,intensity = sum(intensity)),
                                       area          = ddply(peaklist_in,c("Acquisition","Peaklist"),summarise,area = sum(area)),
                                       count         = ddply(peaklist_in,c("Acquisition","Peaklist"),summarise,count = sum(count))
                          ),
                          max = switch(plot_var,
                                       m.z           = ddply(peaklist_in,c("Acquisition","Peaklist"),summarise,m.z = max(m.z)),
                                       SN            = ddply(peaklist_in,c("Acquisition","Peaklist"),summarise,SN = max(SN)),
                                       QualityFactor = ddply(peaklist_in,c("Acquisition","Peaklist"),summarise,QualityFactor = max(QualityFactor)),
                                       Resolution    = ddply(peaklist_in,c("Acquisition","Peaklist"),summarise,Resolution = max(Resolution)),
                                       intensity     = ddply(peaklist_in,c("Acquisition","Peaklist"),summarise,intensity = max(intensity)),
                                       area          = ddply(peaklist_in,c("Acquisition","Peaklist"),summarise,area = max(area))
                          )
    )
  }

  if (is.na(match(plot_var_transform,c("none","log")))){
    print(paste("Transformation",plot_var_transform,"is not currently supported."))
    print("Defaulting to none")
    plot_var_transform = "none"
    print("In order to use non-standard transformations the spatialPlot() function will need to be modified.")
  }
  if (plot_var_transform == "log"){
    peaklist_in[,plot_var] = round(log(1 + peaklist_in[,plot_var]),5)
  }


  p <- switch(plot_var,
              m.z           = ggplot(peaklist_in,aes(x=Acquisition,y=m.z)),
              SN            = ggplot(peaklist_in,aes(x=Acquisition,y=SN)),
              QualityFactor = ggplot(peaklist_in,aes(x=Acquisition,y=QualityFactor)),
              Resolution    = ggplot(peaklist_in,aes(x=Acquisition,y=Resolution)),
              intensity     = ggplot(peaklist_in,aes(x=Acquisition,y=intensity)),
              area          = ggplot(peaklist_in,aes(x=Acquisition,y=area)),
              count         = ggplot(peaklist_in,aes(x=Acquisition,y=count))
  )
  p <- p + layer(geom = "point",alpha=I(1/12))
  p <- p + layer(geom="smooth",method="gam",formula=y~s(x, bs="cs"),size=I(2),level=0.99)
  p <- p + ggtitle(plot_name)


  if (plot_name != ""){
    plot_name = paste(plot_name,"_",sep="")
  }

  if (save_plot){
    ggsave(paste(plot_name,paste(mult_peaks,plot_var,"transformation",plot_var_transform,sep="_"),".png",sep=""),p)
  }

  return(p)

}





