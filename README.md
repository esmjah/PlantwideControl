# Introduction 
This repository contains the codes for simulation is article "Plantwide control of an oil production network"m submitten to [Computers & Chemical Engineering](https://www.journals.elsevier.com/computers-and-chemical-engineering). The article is under under review. A draft version of the article can be seen [here](Article/PlanWideOil_June2019.pdf).

# Models
### Simplified dynamic model:
The simplified dynamic model of the process is implemented in Modelica. The model components have been developed and testes using the [Dymola](https://www.3ds.com/products-services/catia/products/dymola/) software. The models used in our simulation are availabe the [networkModels](networkModels/) directory.

The Modelica model consists of two main components:


# CasADi
The Modelica compiler generates a functional mock-up unit (FMU), which is a standard model component that can be shared with other applications. The resulting model was imported to CasADi
[(Andersson et al., 2019)](http://www.optimization-online.org/DB_FILE/2018/01/6420.pdf) 
which includes efficient automatic
differentiation techniques. NMPC and EKF were imple-
mented using the CasADi verion 2.0.0. The communica-
tion between the Olga Simulator and the controllers was
500 done by OPC Data Access where the Olga OPC Server is
a built-in module of the simulator and the OPC client is
coded in Python.
