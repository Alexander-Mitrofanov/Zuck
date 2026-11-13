import {cp,mkdir,readFile,writeFile} from 'node:fs/promises';
import {makeApi} from '../api.js';
import {createHash} from 'node:crypto';
const hash = source => createHash('sha256').update(source).digest('hex').slice(0,12);
const base=process.env.ZUCK_API_BASE_URL;
if (!base || !base.startsWith('https://')) throw new Error('ZUCK_API_BASE_URL must be an exact HTTPS API address.');
makeApi(base);
await mkdir('dist',{recursive:true});
for (const file of ['index.html','style.css','app.js','api.js','assets','examples','LICENSE']) await cp(file,`dist/${file}`,{recursive:true});
await writeFile('dist/config.js',`window.ZUCK_CONFIG = ${JSON.stringify({apiBase:base})};\n`);
await writeFile('dist/.nojekyll','');
const apiHash = hash(await readFile('dist/api.js'));
const app = (await readFile('dist/app.js','utf8')).replace('./api.js', `./api.js?v=${apiHash}`);
await writeFile('dist/app.js', app);
let index = (await readFile('dist/index.html','utf8')).replace(' http://127.0.0.1:8001','');
for (const file of ['style.css','config.js','app.js']) {
  const version = hash(await readFile(`dist/${file}`));
  index = index.replace(`"${file}"`, `"${file}?v=${version}"`);
}
await writeFile('dist/index.html', index);
console.log('Static Zuck frontend ready.');
