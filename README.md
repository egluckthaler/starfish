# starfish :rocket:

<img
  src="/assets/element_logo_color.png"
  align = "right"
  style="margin: 0 auto; width: 182px; height: 200px">

```starfish``` is a modular toolkit for large mobile element discovery. Built primarily for annotating [giant *Starship* elements](https://academic.oup.com/mbe/article/39/5/msac109/6588634) in fungal genomes, it can be easily adapted to find any large mobile element (≥6kb) that shares the same basic architecture as a fungal *Starship* or a bacterial integrative and conjugative element: a "captain" gene with zero or more "cargo" genes downstream of its 3' end.

## Overview

The ```starfish``` workflow is organized into three main modules: Gene Finder, Element Finder, and Region Finder. Each has dedicated commands that are typically run sequentially. Auxiliary commands that provide additional utilities and generate visualizations are also available through the commandline. Several useful stand-alone scripts are located in the `/scripts` directory. 

<img
  src="/assets/starfishWorkflow.png"
  style="display: center; margin: 0 auto; max-width: 400px">

## Documentation

Head to our [GitHub Wiki](https://github.com/egluckthaler/starfish/wiki) for useful resources, including [installation instructions](https://github.com/egluckthaler/starfish/wiki/Installation), a [manual](https://github.com/egluckthaler/starfish/wiki/Manual) with important details and considerations for each command, and a [step-by-step tutorial](https://github.com/egluckthaler/starfish/wiki/Step-by-step-tutorial). If you run into difficulties, please open an issue on [GitHub](https://github.com/egluckthaler/starfish/issues)

## Citations and dependencies

Many ```starfish``` commands have dependencies that are stand-alone programs in their own right. If you use ```starfish``` in your research, please contact us as it has not yet been published.

Upon publication, please cite both the main ```starfish``` manuscript in addition to any dependencies you may have used (see Table below for a guide). For example:
> We used starfish v1.0 in conjunction with metaeuk, mummer4, and blastn to annotate and visualize Starships in our genome assemblies.

| Command | Dependency | Citation |
|:---:|:---:|:---|
|`annotate`, `augment`| `metaeuk`, `hmmer`, `bedtools` | TBD |
|`insert`, `extend`| `blastn`,`mummer4` | TBD |
|`flank`| `cnef` | TBD |
|`sim`| `sourmash` | TBD |
|`group`| `mcl` | TBD |
|`*-viz`|`circos`,`gggenomes`,`mummer4`, <br />`mafft`, `minimap2`| TBD |

## License

```starfish``` is under active development and will soon be published. If you want to use ```starfish``` in your research, please contact us.

~~```starfish``` is an open source tool available under the GNU Affero General Public License version 3.0 or greater.~~
