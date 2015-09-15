CC= g++
CFLAGS=  -Wno-write-strings

MAINS_DIR= Code/
ATMO_DIR=  Code/Atmo/
GEOAC_DIR= Code/GeoAc/

INSTALL_DIR= /usr/local/bin

GeoAc2D:
	${CC}  ${ATMO_DIR}G2S_Spline1D.cpp ${ATMO_DIR}Atmo_State.Absorption.cpp ${GEOAC_DIR}GeoAc.Parameters.cpp ${GEOAC_DIR}GeoAc.EquationSets.2DStratified.cpp ${GEOAC_DIR}GeoAc.Solver.cpp ${GEOAC_DIR}GeoAc.Interface.cpp ${MAINS_DIR}GeoAc2D_main.cpp ${CFLAGS} -o GeoAc2D

GeoAc3D:
	${CC}  ${ATMO_DIR}G2S_Spline1D.cpp ${ATMO_DIR}Atmo_State.Absorption.cpp ${GEOAC_DIR}GeoAc.Parameters.cpp ${GEOAC_DIR}GeoAc.Eigenray.cpp ${GEOAC_DIR}GeoAc.EquationSets.3DStratified.cpp ${GEOAC_DIR}GeoAc.Solver.cpp ${GEOAC_DIR}GeoAc.Interface.cpp ${MAINS_DIR}GeoAc3D_main.cpp ${CFLAGS} -o GeoAc3D 

GeoAc3D.RngDep:
	${CC}  ${ATMO_DIR}G2S_MultiDimSpline3D.cpp ${ATMO_DIR}Atmo_State.Absorption.cpp ${GEOAC_DIR}GeoAc.Parameters.RngDep.cpp ${GEOAC_DIR}GeoAc.Eigenray.cpp ${GEOAC_DIR}GeoAc.EquationSets.3DRngDep.cpp ${GEOAC_DIR}GeoAc.Solver.cpp ${GEOAC_DIR}GeoAc.Interface.cpp ${MAINS_DIR}GeoAc3D.RngDep_main.cpp ${CFLAGS} -o GeoAc3D.RngDep 

GeoAcGlobal:
	${CC}  ${ATMO_DIR}G2S_GlobalSpline1D.cpp ${ATMO_DIR}Atmo_State.Absorption.Global.cpp ${GEOAC_DIR}GeoAc.Parameters.Global.cpp ${GEOAC_DIR}GeoAc.Eigenray.Global.cpp ${GEOAC_DIR}GeoAc.EquationSets.Global.cpp ${GEOAC_DIR}GeoAc.Solver.cpp ${GEOAC_DIR}GeoAc.Interface.Global.cpp ${MAINS_DIR}GeoAcGlobal_main.cpp ${CFLAGS} -o GeoAcGlobal

GeoAcGlobal.RngDep:
	${CC}  ${ATMO_DIR}G2S_GlobalMultiDimSpline3D.cpp ${ATMO_DIR}Atmo_State.Absorption.Global.cpp ${GEOAC_DIR}GeoAc.Parameters.Global.cpp ${GEOAC_DIR}GeoAc.Eigenray.Global.cpp ${GEOAC_DIR}GeoAc.EquationSets.GlobalRngDep.cpp ${GEOAC_DIR}GeoAc.Solver.cpp ${GEOAC_DIR}GeoAc.Interface.Global.cpp ${MAINS_DIR}GeoAcGlobal.RngDep_main.cpp ${CFLAGS} -o GeoAcGlobal.RngDep

all: GeoAc2D GeoAc3D GeoAc3D.RngDep GeoAcGlobal GeoAcGlobal.RngDep

clean: 
	rm GeoAc2D GeoAc3D GeoAc3D.RngDep GeoAcGlobal GeoAcGlobal.RngDep

install: GeoAc2D GeoAc3D GeoAc3D.RngDep GeoAcGlobal GeoAcGlobal.RngDep
	install -m 0755 GeoAc2D ${INSTALL_DIR}
	install -m 0755 GeoAc3D ${INSTALL_DIR}
	install -m 0755 GeoAc3D.RngDep ${INSTALL_DIR}
	install -m 0755 GeoAcGlobal ${INSTALL_DIR}
	install -m 0755 GeoAcGlobal.RngDep ${INSTALL_DIR}

uninstall:
	rm ${INSTALL_DIR}/GeoAc2D 
	rm ${INSTALL_DIR}/GeoAc3D 
	rm ${INSTALL_DIR}/GeoAc3D.RngDep 
	rm ${INSTALL_DIR}/GeoAcGlobal 
	rm ${INSTALL_DIR}/GeoAcGlobal.RngDep 
	rm GeoAc2D GeoAc3D GeoAc3D.RngDep GeoAcGlobal GeoAcGlobal.RngDep