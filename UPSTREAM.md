# Authors and dependencies

Zuck includes detection code by Fikrat Talibli and Björn Voß, with local
modifications made in November 2026. The source is distributed under the
GNU General Public License version 3; see `LICENSE.txt`. Existing copyright
and license notices in bundled dependencies remain in their source trees.

The development starting revision of the detection core was
`26a6246087b56ad14e7cf01e9b8d2fdb0217484c`. Bundled dependency revisions are:

| Component | Revision |
|---|---|
| cft | `6cb63d0036223b8be9c4d01f95b296d641209f95` |
| kseqpp | `50b7fecb7dbd10758600a6f5d4fef6b2f3972f79` |
| megahit | `5f329c6d24a1480b75145a4c14567a25453b95bf` |
| spoa | `167bdcc41a25e397cbd4487b6338989266d815d8` |

Optional vendored test/data collections are not included. CMake also fetches
parallel-hashmap v1.3.11, fmt 10.2.1, and zlib v1.3.1 when needed.

The project rename does not claim that existing code was independently
authored. Original documentation and development records are preserved in
the workspace's frozen provenance archive.
