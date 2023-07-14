using LocalRegistry
using Pkg
using SorptionModels


Pkg.Registry.update()
register(SorptionModels; registry="MembraneRegistry")