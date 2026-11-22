"""Recover tied observed terminal alleles without replacing baseline evidence."""
from collections import Counter
import csv,json,random,subprocess,sys,tempfile
from pathlib import Path


def read_rows(folder):
    with (folder/'entries.tsv').open() as handle:
        return Counter((r['repeat_variant'],r['spacer']) for r in csv.DictReader(handle,delimiter='\t'))


def main():
    binary=Path(sys.argv[1]).resolve()
    rng=random.Random(872381)
    dna='ACGT'
    sequence=lambda n:''.join(rng.choice(dna) for _ in range(n))
    repeat=sequence(30)
    # Balanced spacer ends prevent accidental repeat-boundary extension. Truth
    # consists of complete generating pairs, independent of inferred consensus.
    independent=[dna[i%4]+sequence(30)+dna[(i//4+i)%4] for i in range(12)]
    anchor=independent[0]
    shorter=anchor[:-1]
    expected_new=(repeat,shorter)
    with tempfile.TemporaryDirectory(prefix='zuck-terminal-rescue-') as tmp:
        root=Path(tmp);serial=0
        def run(anchor_support,short_support,rescue,mode='two-sided',settings=False):
            nonlocal serial
            serial+=1
            rows=[repeat+s for s in independent[1:]]+[repeat+anchor]*anchor_support+[repeat+shorter]*short_support
            cycles=root/f'cycles-{serial}.txt';cycles.write_text('\n'.join(c+c[:22] for c in rows)+'\n')
            out=root/f'out-{serial}'
            command=[str(binary),'--cycles',str(cycles),'--output-folder',str(out),'--threads','1','--ram','1G','--boundary-mode',mode]
            if settings:
                config=root/f'settings-{serial}.txt';config.write_text('terminal-variant-rescue=true\n')
                command+=['--settings',str(config)]
            if rescue is not None:command+=['--terminal-variant-rescue',str(rescue).lower()]
            result=subprocess.run(command,capture_output=True,text=True,timeout=30)
            if result.returncode:raise AssertionError(result.stdout+result.stderr)
            return read_rows(out),json.loads((out/'run_manifest.json').read_text())
        baseline,off=run(3,3,False)
        recovered,on=run(3,3,True)
        assert expected_new not in baseline, 'fixture must expose representative suppression'
        assert recovered-baseline==Counter({expected_new:1})
        assert not baseline-recovered, 'rescue removed/replaced an original oriented pair'
        assert on['counts']['terminal_variant_rows_rescued']==1
        assert off['counts']['terminal_variant_rows_rescued']==0
        assert on['settings']['terminal_variant_rescue'] is True
        assert off['settings']['terminal_variant_rescue'] is False
        for support in [(3,2),(2,2)]:
            before,_=run(*support,False);after,manifest=run(*support,True)
            assert before==after
            assert manifest['counts']['terminal_variant_rows_rescued']==0
        before,_=run(3,3,False,mode='legacy');after,_=run(3,3,True,mode='legacy')
        assert before==after, 'legacy inference must remain unchanged'
        configured,_=run(3,3,None,settings=True)
        overridden,_=run(3,3,False,settings=True)
        assert configured==recovered and overridden==baseline
        print('Terminal rescue preserves baseline pairs and adds observed tied allele')


if __name__=='__main__':main()
