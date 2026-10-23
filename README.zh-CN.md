# Nivat 星形构形的精确复杂度渐近与有限窗口下界

指定比赛论文中的全部星形构形，在任意固定、有内部的紧凸窗口形状上，都有精确的复杂度二次首项。系数是方向带两两交的并集面积，证明同时给出匹配上界，覆盖 Case A 与 Case B。

## 对象与结论

对象是 Apex Intelligence 于 2026 年 9 月 12 日发布的比赛版 *The Convex Nivat Conjecture: A Complexity Lower Bound for Star Configurations, and a Reduction from Low Convex Complexity to Star Configurations*。PDF SHA-256：`7fd67831155f4226c010fd6af32de558b76643771e6be12dea4eebafa21745a8`。引用采用印刷页码与编号。本包支持针对星形构形的突破贡献，类别为 `supplement`。

设 $`B_i(P)`$ 是固定凸体 $`P`$ 沿分量方向 $`v_i`$ 的两条支撑线之间的完整带，令

```math
\mathcal M(P)=\bigcup_{i\lt j}(B_i(P)\cap B_j(P)).
```

- **定理 1：**全部 CN 星形都满足 $`P_\theta(nP\cap\mathbb Z^2)=\mathrm{area}(\mathcal M(P))n^2+O_{\theta,P}(n)`$。$`P`$ 可为任意固定、有内部的紧凸集。
- **定理 2：**全部星形、任意有限非空窗口的显式下界，计入完整标记可见的中心锚点及两方向可恢复的远端锚点。证书常数有有限构造，但小窗口上可能弱于 CN 的 $`|S|+1`$。
- **定理 3：**三线二值例在标准三角窗口上满足 $`P_c(S_n)=31n^2/4+11n/4+O(1)`$。这里只确定该例的一次项。

同一 496 点窗口 $`S_{30}`$ 上，三线例下界为 4059，CN 给出 497，真实模式数为 7037；Case A 象限例下界为 786，CN 同为 497，真实模式数为 1022。

## 阅读顺序

先读 [对象与定理](proof/objects.md)，再按以下顺序阅读英文完整证明：

1. [条带信号与锚点恢复](proof/signals-and-anchors.md)：Case B 的非零条带和极端行恢复。
2. [首项几何与匹配上界](proof/leading-term.md)：固定侵蚀损失、格点计数、全部星形的上界及有限容斥面积计算。
3. [全部星形的下界](proof/all-star-extension.md)：两端共同极值、有界核心、远端恢复及完整标记观察。
4. [显式有限窗口界](proof/finite-window-bound.md)：定理 2、自动常数构造和象限例。
5. [三线例的一次项](proof/three-line-second-order.md)：定理 3 与碰撞分类。
6. [比较、归因与范围](proof/remarks-and-scope.md)：4059 的计算、扫描矩形和象限阈值单元的覆盖证明，以及经典尖锐性备注。

[文献表](sources/bibliography.md)列明实际核对的版本与页码。差分观测中心属于窗口凸包的依据，以及原点系数抵消的反例，都已明确写入证明。

## 复现

只需 Python 3.9 及以上版本的标准库。在仓库根目录运行：

```sh
python3 -B scripts/verify.py
```

程序依次重跑四份核验，将生成的 JSON 与随包确定性输出逐项比较；每份程序只写自身旁边的 JSON，不请求网络或读取账户配置。四份程序分别检查精确面积与锚点数、一般星形两端信号、所有 $`3\le n\le164`$ 的三线模式及碰撞、有限窗口与证书压力样例。[三线扫描矩形](proof/remarks-and-scope.md#11-complete-line-enumeration)与[象限全部阈值单元](proof/remarks-and-scope.md#12-complete-quadrant-enumeration)的覆盖论证均在证明内；程序完整目录见[英文 README](README.md#reproducing-the-computations)。有限计算不替代一般证明。

可用 `shasum -a 256 -c MANIFEST.sha256` 核对整个包的指纹。

## 归因与限制

CN 提供星形定义、差分几何、非零条带论证及平移条带乘积标记。这里用极端行、共同端部极值与完整标记观察证明全局唯一性，再结合上下界得到精确系数。证明不使用 Kari–Szabados 理想结构、谱多项式或标量理想桥。

二次增长本身并非新贡献：Kari–Szabados 已有矩形窗口的面积阶下界，部分证明还给出二次增量。本包的定量结论是全部星形、任意固定凸形状的精确系数及匹配上界。未做全面文献检索；单点缺陷只列为经典尖锐性备注，历史归因经 Szabados 引言及其所引 Cassaigne 文献说明。

不声称提高任意构形的周期性阈值，也不声称一般二阶公式、最终准多项式或证明助手认证。

## AI 声明

GPT-6.1 Sol, in OpenAI Codex under human direction, developed the proof, ran the computations, and drafted this text; separate GPT-6.1 Sol sessions in fresh contexts performed adversarial reviews. Claude planned the work, checked key steps against the sources, and reviewed and edited the final text. GPT-6 Astra, also in OpenAI Codex, advised on research direction and submission scope. No human expert has certified the work.

## 许可

原创散文及编排采用 CC BY 4.0；`scripts/` 代码采用 MIT。外部论文保留其原有权利，仅引用而不随包分发。详见 [LICENSE](LICENSE)。
