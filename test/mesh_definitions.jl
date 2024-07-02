const global Cell_ = Dimension{:Cell_, HORIZONTAL}
const global K_ = Dimension{:K_, HORIZONTAL}
const global Edge_ = Dimension{:Edge_, HORIZONTAL}
const global Vertex_ = Dimension{:Vertex_, HORIZONTAL}
const global V2VDim_ = Dimension{:V2VDim_, LOCAL}
const global V2EDim_ = Dimension{:V2EDim_, LOCAL}
const global E2VDim_ = Dimension{:E2VDim_, LOCAL}
const global E2CDim_ = Dimension{:E2CDim_, LOCAL}
const global C2EDim_ = Dimension{:C2EDim_, LOCAL}
const global Cell = Cell_()
const global K = K_()
const global Edge = Edge_()
const global Vertex = Vertex_()
const global V2VDim = V2VDim_()
const global V2EDim = V2EDim_()
const global E2VDim = E2VDim_()
const global E2CDim = E2CDim_()
const global C2EDim = C2EDim_()

const global V2V = FieldOffset("V2V", source=Vertex, target=(Vertex, V2VDim))
const global E2V = FieldOffset("E2V", source=Vertex, target=(Edge, E2VDim))
const global V2E = FieldOffset("V2E", source=Edge, target=(Vertex, V2EDim))
const global E2C = FieldOffset("E2C", source=Cell, target=(Edge, E2CDim))
const global C2E = FieldOffset("C2E", source=Edge, target=(Cell, C2EDim))
const global Koff = FieldOffset("Koff", source=K, target=K)