###
# NADA for R by Lopaka Lee.
#
# Version 1.2-1
# Copyright (2004, 2005) Lopaka Lee
#
# A S-language software module based on 
# methodologies described by Dennis R. Helsel in his book 
# Nondetects and Data Analysis: Statistics for Censored Environmental Data.
#
# NADA is free software; you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free
# Software Foundation; either version 2 of the License, or (at your option)
# any later version.
# 
# NADA is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
# for more details.  You should have received a copy of the GNU General
# Public License along with NADA; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
###

### Definitions common to All sections in this package

## Globals

.onLoad = function(lib, pkg) require(methods)

NADAprobs = c(0.05,0.10,0.25,0.50,0.75,0.90,0.95)

## Generics

setGeneric("plot", function(x, y, ...) standardGeneric("plot"))

setGeneric("print", function(x, ...) standardGeneric("print"))

setGeneric("summary", function(object, ...) standardGeneric("summary"))

setGeneric("mean", function(x, ...) standardGeneric("mean"))

setGeneric("sd", function(x, na.rm=FALSE) standardGeneric("sd"))

setGeneric("median", function(x, na.rm=FALSE) standardGeneric("median"))

setGeneric("quantile", function(x, ...) standardGeneric("quantile"))

setGeneric("predict", function(object, ...) standardGeneric("predict"))

setGeneric("lines", function(x, ...) standardGeneric("lines"))

## Classes

setClass("NADAlist", "list")

## Methods

setMethod("print", signature("NADAlist"), function(x, ...)
{
    tag = names(x)
    for (i in 1:length(x))
      {
        cat(tag[i], "\n")
        print(x[[i]])
        cat("\n")
      }
})

#-->> BEGIN general utility functions

##
# splitQual extracts qualifed and unqualifed vectors from a vector
# containing concatenated qualifiying characters and numeric values
# like "<0.5".  Only handles one kind of censoring character/symbol.
splitQual =
function(v, qual.symbol = "<")
{
    qual.index = grep(qual.symbol, x=as.character(v))  

    qual.chars = as.character(v[qual.index])
    qual = as.numeric(sub(qual.symbol, "", qual.chars))

    unqual.index = -1 * qual.index
    unqual = as.numeric(as.character(v[unqual.index]))
    
    return(list(qual         = qual, 
                unqual       = unqual, 
                qual.index   = qual.index,
                unqual.index = unqual.index))
}

## pct.censored  -- Simple function to save some typing
pct.censored =
function(obs, censored)
{
    if (!is.logical(censored)) 
      {
        stop("censored indicator must be logical vector!\n")
      }

    return(100*(length(obs[censored])/length(obs)))
}

#-->> END general utility functions

