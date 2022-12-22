# starfish: the STARship FInder SHell :rocket:

starfish is a computational workflow for large mobile element discovery. Built primarily for annotating [giant *Starship* elements](https://academic.oup.com/mbe/article/39/5/msac109/6588634) in fungal genomes, it can be easily adapted to find any large mobile element (≥6kb) that shares the same basic architecture as a fungal *Starship* or bacterial integrative and conjugative element: a "captain" gene with zero or more "cargo" genes located downstream of its 3' end.

## overview

starfish is organized into three main modules: gene finder, element finder, and region finder. Each has several dedicated commands that are typically run sequentially. Many auxiliary commands that provide additional analyses and generate visualizations are also available through the commandline. Several useful stand-alone scripts are located in the `/scripts` directory. 

<img
  src="/data/starfishWorkflow.jpg"
  alt="starfish workflow and commands"
  title="starfish workflow and commands"
  style="display: center; margin: 0 auto; max-width: 300px">

## installation

### main workflow and dependencies

start by cloning the latest version of this github repository:

```
git clone https://github.com/egluckthaler/starfish.git
cd starfish/
```

then use conda to install other dependencies into a new, dedicated environment:

```
conda install --file meta.yaml
```

activate the new conda environment anytime you want to use starfish:

```
conda activate starfish
```

finally, add the following lines to your $HOME/.bash_profile:

```
export PATH=$PATH:$HOME/starfish/
export PATH=$PATH:$HOME/starfish/CNEFinder/
```

starfish comes with a command auto-completion, which can be activated by adding these additional lines to your $HOME/.bash_profile:

```
<pre>
if [ -f $HOME/starfish/util/bash-completion.sh ]; then
    source $HOME/starfish/util/bash-completion.sh
fi
</pre>
```

and you're done!

### additional dependencies

The R package `gggenomes` is used in some visualization commands and can be installed by visiting [the gggenomes github](https://github.com/thackl/gggenomes)

## documentation

The starfish user manual is available through our [GitHub Wiki](https://github.com/egluckthaler/starfish/wiki). The wiki also contains step-by-step [tutorials](https://github.com/egluckthaler/starfish/wiki/Tutorials) to learn how to use starfish with real data. If you run into difficulties, please open an issue on [GitHub](https://github.com/egluckthaler/starfish/issues)

## citations and dependencies

Some starfish commands have dependencies that are stand-alone programs in their own right. If you use starfish in your research, please cite both our forthcoming publication and any  dependencies you may have used (see Table below for a guide). For example:
> We used starfish v1.0 in conjunction with metaeuk, mummer4, and blastn to annotate and visualize Starships in our genome assemblies.

| Command | Dependency | Citation |
|:---:|:---:|:---|
|`annotate`, `augment`| `metaeuk`, `hmmer`, `bedtools` | TBD |
|`insert`, `extend`| `blastn`,`mummer4` | TBD |
|`flank`| `cnef` | TBD |
|`sim`| `sourmash` | TBD |
|`group`| `mcl` | TBD |
|`*-viz`|`circos`,`gggenomes`,`mummer4`, <br />`mafft`, `minimap2`| TBD |

## license

starfish is an open source tool available under the GNU Affero General Public License version 3.0 or greater.