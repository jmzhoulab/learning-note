# LOUVAIN算法模型

​		Louvain算法是一种基于多层次（逐轮启发式迭代）优化Modularity的算法。Modularity函数最初被用于衡量社区发现算法结果的质量，它能够刻画发现的社区的紧密程度。

同时，Modularity函数既然能刻画社区的紧密程度，也就能够被用来当作一个优化函数（目标函数），即将结点加入它的某个邻居所在的社区中，如果能够提升当前社区结构的modularity。则说明这次迭代优化是可接受的。

下面我们来讨论Louvain算法模型的核心组件。

## Modularity的定义

​		模块度是评估一个社区网络划分好坏的度量方法，它的物理含义是社区内节点的连边数与随机情况下的边数只差，它的取值范围是 [−1/2,1)，其定义如下：
$$
Q=\frac{1}{2m}\sum_{i,j}\bigg[A_{ij}-\frac{k_ik_j}{2m}\bigg]\delta(c_i, c_j)
$$
其中，$A_{ij}$是节点$i$到节点$j$关系边的权重，$k_i=\sum_{j}A_{ij}$是链接节点$i$的边的权重和，$c_i$是节点$i$所在社区的标签，$\delta(u,v)$是示性函数，即，当$u=v$时，$\delta(u,v)=1$，否则$\delta(u,v)=0$。$m=\frac{1}{2}\sum_{ij}A_{ij}$

​		模块度的公式可以进行如下简化：
$$
\begin{align}
\nonumber Q &= \frac{1}{2m}\sum_{i,j}\bigg[A_{ij}-\frac{k_ik_j}{2m} \bigg]\delta(c_i, c_j) \\
\nonumber  &= \frac{1}{2m}\bigg[\sum_{i,j}A_{ij}-\frac{\sum_{i}k_i\sum_{j}k_j}{2m}\bigg] \delta(c_i, c_j)\\
  &=\frac{1}{2m}\sum_{c}\bigg[\Sigma_{in}-\frac{(\Sigma_{tot})^2}{2m}\bigg]
\end{align}
$$
其中 $\Sigma_{in}$ 表示社区 $c$内的边的权重之和； $\Sigma_{tot}$表示与社区$c$内的节点相连的所有边的权重之和。

​		上面的公式还可以进一步简化成:
$$
\begin{align}
\nonumber Q &=\frac{1}{2m}\sum_{c}\bigg[\Sigma_{in}-\frac{(\Sigma_{tot})^2}{2m}\bigg] \\
  &=\sum_{c}[e_{c}-a_{c}^2]
\end{align}
$$
**模块度的直观理解：**

​		首先modularity是针对一个社区的所有节点进行了累加计算。其计算公式背后体现了这种思想：社区内部边的权重减去所有与社区节点相连的边的权重和，对无向图更好理解，即社区内部边的度数减去社区内节点的总度数。

​		**可以直观去想象一下，如果一个社区节点完全是“封闭的（即所有节点都互相内部连接，但是不和社区外部其他节点有连接，则modularity公式的计算结果为1）**



## 模块度增量 $\Delta Q$

​		模块增益度是评价本次迭代效果好坏的数值化指标，这是一种启发式的优化过程。**类似决策树中的熵增益启发式评价**。

​		当一个孤立节点i（即节点i还未被分配到任何社区）加入到社区$C$的模块度增益计算公式如下:
$$
\begin{align}
\nonumber \Delta Q &=\bigg[\frac{\Sigma_{in}+k_{i,in}}{2m} - \bigg(\frac{\Sigma_{tot}+k_i}{2m}\bigg)^2\bigg] - \bigg[\frac{\Sigma_{in}}{2m} - \bigg(\frac{\Sigma_{tot}}{2m}\bigg)^2
- \bigg(\frac{k_i}{2m}\bigg)^2
\bigg] \\
&=\frac{1}{2m}\bigg(k_{i,in} - \frac{k_i\Sigma_{tot}}{m}\bigg)
\end{align}
$$
其中， $\Sigma_{in}$ 表示社区 $C$内的边的权重之和；$k_{i,in}$时节点i与社区C内节点的权重和； $\Sigma_{tot}$表示与社区$C$内的节点相连的所有边的权重之和。$k_i=\sum_{j}A_{ij}$是链接节点$i$的边的权重和

​		当社区$C_{1}$中的节点$i$重新划分到社区$C_2$时，模块度增益计算公式如下：
$$
\begin{align}
\nonumber \Delta Q =& \frac{1}{2m}\bigg[\Sigma^{(2)}_{in}+k^{(2)}_{i,in} - \frac{(\Sigma^{(2)}_{tot}+k_i)^2}{2m} + \Sigma^{(1)}_{in}-k^{(1)}_{i,in} - \frac{(\Sigma^{(1)}_{tot}-k_i)^2}{2m} \bigg] \\
\nonumber &- \frac{1}{2m}\bigg[\Sigma^{(2)}_{in} - \frac{(\Sigma^{(2)}_{tot})^2}{2m}
+ \Sigma^{(1)}_{in} - \frac{(\Sigma^{(1)}_{tot})^2}{2m}
\bigg] \\
=&\frac{1}{2m}\bigg[k^{(2)}_{i,in} - k^{(1)}_{i,in} - \frac{k_i(\Sigma^{(2)}_{tot}-\Sigma^{(1)}_{tot})}{m}\bigg]
\end{align}
$$
其中， $\Sigma^{l}_{in}$ 表示社区 $C_{l}$内的边的权重之和；$k^l_{i,in}$是节点i与社区$C_l$内节点的权重和； $\Sigma^l_{tot}$表示与社区$C$内的节点相连的所有边的权重之和。$k_i=\sum_{j}A_{ij}$是链接节点$i$的边的权重和

