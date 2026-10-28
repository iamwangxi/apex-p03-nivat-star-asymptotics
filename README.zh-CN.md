# Nivat 星形的精确复杂度渐近与最终拟多项式

全部星形在任意固定、有内部的紧凸形状上都有精确的面积首项；有理多边形上最终恰为二次拟多项式；满支撑非零周期线类还有直接的二阶几何公式。

## 已复核结果与一般系数补充

对象为 Apex Intelligence 的 2026 年 9 月 12 日比赛版 *The Convex Nivat Conjecture*。源文件指纹与复核范围见[状态记录](sources/review-status.md)，文献版本、定理和页码见[文献表](sources/bibliography.md)。

- 全部星形、任意固定凸体：$`P_\theta(nP\cap\mathbb Z^2)=\Lambda(P)n^2+O(n)`$，其中 $`\Lambda(P)`$ 是方向带两两交之并的面积；另有显式有限窗口下界。
- 全部星形、有理多边形：最终精确为 $`\Lambda(P)n^2+\beta_{n\bmod M}n+\gamma_{n\bmod M}`$。这是已有 Presburger 计数理论的星形应用，具体增量在于全部星形的有限描述与面积项识别。
- 满支撑周期线类：非零周期词、零背景、任意有理多边形，包括相交和颜色抵消；直接公式为 $`\beta_r=\mathfrak B_r-\Gamma_r`$。整格多边形的一次项恒定；窗口顶点分母模数只保证该一次项周期，不冒充完整余项的周期。
- [一般系数与充分周期补充页](proof/explicit-coefficients.md)：**经一路全新上下文对抗复核成立**。经过有限消元，以代表多边形的面积一次项、格边界和一维长度求和。它是算法式几何公式，可能非常大，未得到原始方向带上的简短闭式。其余结果均不依赖此页。
- [精确常数项与最小周期补充页](proof/minimal-periods-and-constant-term.md)：**经一路全新上下文对抗复核成立**（2026-10-26）。由顶点速度分母得到更小的充分周期 $`M=qE`$；给出常数项 $`\gamma`$ 的显式有限求和；用分圆多项式整除精确判定 $`\beta`$、$`\gamma`$ 与整个拟多项式的最小周期；三线着色族满足 $`\beta=3-\tfrac14\mathbf 1[c_0=c_2]`$，排除只用无色几何或分量独立相加的公式。

本次整合与一般系数补充已于 2026-10-25、常数项与最小周期补充已于 2026-10-26 各经一路全新上下文 GPT-6.1 Sol 对抗复核，未发现致命或证明修补问题；这是 AI 复核，无人类专家或证明助手认证。指纹与范围见[状态记录](sources/review-status.md)。

## 完整证明

从[对象定义](proof/objects.md)开始，依次可读[信号与锚点](proof/signals-and-anchors.md)、[首项几何](proof/leading-term.md)、[全部星形下界](proof/all-star-extension.md)、[有限窗口界](proof/finite-window-bound.md)、[最终拟多项式](proof/eventual-quasipolynomial.md)、[满支撑二阶公式](proof/full-support-lines.md)。

[三线例](proof/three-line-second-order.md)另给一次项的完整碰撞证明；[数值与边界](proof/remarks-and-scope.md)含扫描矩形、象限阈值单元的全局覆盖证明。满支撑页面保留自己的局部定理和公式编号。

同一三线三角例在 $`n=30`$ 时有496个站点，下界4059，原论文给497，精确模式数7037。其首项和一次项为 $`31/4`$ 与 $`11/4`$，常数余项最终周期；没有把有限数据拟合出的周期表当作定理。不同的 Case A 象限例在同一窗口上，下界786，精确值1022。

## 复现

只需 Python 3.9 及以上版本的标准库；无需下载、账户或外部论文。仓库根目录执行：

~~~sh
python3 -B scripts/verify.py
~~~

程序重跑8份诊断，比较7份JSON和1份文本的完整确定性输出。保留原4份程序，新增满支撑公式、有限描述样例及一般系数补充的几何诊断。每份程序只写自身相邻输出。有限计算不替代无限证明；详见[英文说明](README.md#reproducing-the-computations)。

~~~sh
shasum -a 256 -c MANIFEST.sha256
~~~

可核对整个包的指纹。清单不包含它自身。

## 范围与归因

CN 提供差分、条带信号和平移乘积标记。Kari–Szabados 已有面积阶增长结果，二次增长本身不作新颖性主张。拟多项式所用计数原理来自 Bogart–Goodrick–Woods；补充的一般系数构造同样使用经典消元与格点几何，不声称新的通用计数理论。

不声称提高任意构形的周期性阈值，不声称最小周期、全部星形在整格窗口上一次项恒定、证明助手认证或文献优先权。无理窗口上的满支撑十字反例阻止把周期一次项结论直接推广到所有凸窗口。

## AI 声明

以下沿用来源工作的角色声明；本次整合及两份补充已分别于 2026-10-25、10-26 接受全新上下文 GPT-6.1 Sol 对抗复核；Claude 完成最终编排与编辑审阅。

GPT-6.1 Sol, in OpenAI Codex under human direction, developed the proof, ran the computations, and drafted this text; separate GPT-6.1 Sol sessions in fresh contexts performed adversarial reviews. Claude planned the work, checked key steps against the sources, and reviewed and edited the final text. GPT-6 Astra, also in OpenAI Codex, advised on research direction and submission scope. No human expert has certified the work.

## 许可

原创文字及编排采用 CC BY 4.0，`scripts/` 代码采用 MIT。外部论文仅引用，不随包分发。见[许可](LICENSE)。
