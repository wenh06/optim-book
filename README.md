# optim-book

本科生教材《最优化理论与方法》关于二次规划问题的一节（第七章第3节）

主文件是[`chap7.tex`](chap7.tex)，每一个小节（subsection）都单独写在了各自的文件内，这些文件都放在文件夹[`content`](content)内。

文件[`chap7_in_one.tex`](chap7_in_one.tex)是利用[`bib-lookup`](https://pypi.org/project/bib-lookup/)自动生成的，
具体可见[`adjust-formats.ipynb`](notebooks/adjust-formats.ipynb)，也可直接在项目的根文件夹下执行下述命令

```python
bib-lookup --gather ./chap7.tex
```

文件[`chap7.md`](chap7.md)是利用[`pandoc`](https://pandoc.org/)执行如下命令自动生成的Markdown文件

```bash
pandoc chap7.tex -o chap7.md
```
