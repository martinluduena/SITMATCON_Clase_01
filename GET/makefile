#!/bin/bash


        # COMPILADOR
	FC=gfortran

  	 DEFS =-O0
	 #DEFS +=-fbacktrace
	 #DEFS +=-fdump-parse-tree-fdump-fortran-original
	 #DEFS +=-fdump-core-fdump-fortran-optimize
	 DEFS +=-pedantic
	 DEFS +=-ffree-line-length-none 

	# LIBRERIAS QUE UTILICE EL PROGRAMA
	LIB=-L/usr/local/lib -lfftw3 -lm

        # LISTA DE OBJETOS
 	MODULES = mod_count_lines.o
	MODULES += mod_ran2.o
	MODULES += mod_RDF.o
	MODULES += mod_Clean.o
	MODULES += mod_SMTB.o
	MODULES += biseccion.o
	MODULES += Simpson.o
	MODULES += weighting_functions.o
	MODULES += mod_DoS_2PT.o
	MODULES += mod_low_DoS.o
	MODULES += mod_FFT.o
	MODULES += mod_FEFF.o
	MODULES += mod_BORDER.o
	MODULES += mod_PDF.o
	MODULES += mod_get_comand.o
        
	# PROGRAMA PRINCIPAL
	MAIN = get_command.o                   
	
make:

	@${FC} -c ${DEFSdeb}  mod_count_lines.f90
	@${FC} -c ${DEFSdeb}  mod_ran2.f90
	@${FC} -c ${DEFSdeb}  mod_RDF.f90                                                     
	@${FC} -c ${DEFSdeb}  mod_Clean.f90                                                   
	@${FC} -c ${DEFSdeb}  mod_SMTB.f90                                                    
	@${FC} -c ${DEFSdeb}  biseccion.f90
	@${FC} -c ${DEFSdeb}  Simpson.f90
	@${FC} -c ${DEFSdeb}  weighting_functions.f90
	@${FC} -c ${DEFSdeb}  mod_DoS_2PT.f90   
	@${FC} -c ${DEFSdeb}  mod_low_DoS.f90                                                  
	@${FC} -c ${DEFSdeb}  mod_FFT.f90                               
	@${FC} -c ${DEFSdeb}  mod_FEFF.f90                                                    
	@${FC} -c ${DEFSdeb}  mod_BORDER.f90                                                  
	@${FC} -c ${DEFSdeb}  mod_PDF.f90
	@${FC} -c ${DEFSdeb}  mod_get_comand.f90
	@${FC} -c ${DEFSdeb}  get_command.f90                   
	
	@echo "=================== COMPILACION: DEBUG ============"
	
	${FC} -o get ${MODULES} ${MAIN} ${LIB}
	
	@rm -f *.o
	@rm -f *.mod
