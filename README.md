# starfish: the STARship FInder SHell :rocket:

starfish is a computational workflow for large mobile element discovery. Built primarily for annotating giant [*Starship* elements](https://academic.oup.com/mbe/article/39/5/msac109/6588634) in fungal genomes, it can be easily adapted to find any large mobile element (≥6kb) that shares the same basic architecture as a fungal *Starship* or bacterial integrative and conjugative element: a "captain" gene with zero or more "cargo" genes located downstream of its 3' end.

## overview

starfish is organized into three main modules: gene finder, element finder, and region finder. Each has several dedicated commands that are typically run sequentially. Many auxiliary commands that provide additional analyses and generate visualizations are also available through the commandline. Several useful stand-alone scripts are located in the /scripts/ directory. 

## installation

- the easiest way to install starfish is using git:

```
git clone https://github.com/egluckthaler/starfish.git
cd starfish/
```

- we recommend using conda to install dependencies into a dedicated environment:

```
conda install --file meta.yml/
conda activate starfish
```

- The R package `gggenomes` is used in some visualization commands and can be installed by visiting [the gggenomes github](https://github.com/thackl/gggenomes)

- starfish comes with a command auto-completion, which can be activated by adding the following lines to your $HOME/.bash_profile:

```
<pre>
if [ -f $HOME/scripts/starfish/util/bash-completion.sh ]; then
    source $HOME/scripts/starfish/util/bash-completion.sh
fi
</pre>
```

## documentation

The starfish user manual is available through our [GitHub Wiki](https://github.com/egluckthaler/starfish/wiki). The wiki also contains step-by-step [tutorials](https://github.com/egluckthaler/starfish/wiki/Tutorials) to learn how to use starfish with real data. If you run into difficulties, please open an issue on [GitHub](https://github.com/egluckthaler/starfish/issues)

## citations

Some starfish commands have dependencies that are stand-alone programs in their own right. If you use starfish in your research, please cite both our forthcoming publication and any  dependencies you may have used (see Table below for a guide). For example:
> We used starfish v1.0 in conjunction with metaeuk, mummer4, and blastn to annotate and visualize Starships in our genome assemblies.

| Command | Dependency | Citation |
|:---:|:---:|:---|
|`annotate`, `augment`| `metaeuk`, `hmmer`, `bedtools` | TBD |
|`insert`, `extend`| `blastn`,`mummer4` | TBD |
|`flank`| `cnef` | TBD |
|`*-viz`|`circos`,`gggenomes`,`mummer4`,`mafft`, `minimap2`| TBD |
|`sim`| `sourmash` | TBD |
|`group`| `mcl` | TBD |

## license

