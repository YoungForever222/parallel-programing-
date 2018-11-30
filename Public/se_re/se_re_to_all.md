# 采用MPI_Send和MPI_Recv实现MPI_Allgather的功能
## 算法描述
采用MPI_Send和MPI_Recv实现MPI_Allgather的功能，为了比较其性能，程序以三种方式实现了该功能，传递数据为各个进程所处的节点名称，使用MPI_GET_PROCESSOR_NAME函数获取所在节点名称，然后通过所实现的MPI_Allgather在进程间同步这一数据，并将这一数据用于MPI_Comm_split函数，将进程分组，三种实现如下：
1. 定义主进程，该进程使用MPI_Recv收集其余进程发送的节点名称，组成数组，再将该数组播送给各个进程，实现同步
2. 对于任意进程，无差别地使用MPI_ISEND向全部进程发送自己的节点名称（包括自己），并使用MPI_IRECV收取相同的次数，根据获取数据的节点号填充到数组的相应位置，实现同步
3. 直接调用MPI_Allgather，实现同步
## 程序结构

这里采用Fortran90编写程序,其中：
1. 模块aparallel包含MPI调用过程中需要的参数和变量；
2. 主函数Main包含MPI的初始化，通过参数Use_Se_Re控制采用哪种实现。

## 运行时性能加速