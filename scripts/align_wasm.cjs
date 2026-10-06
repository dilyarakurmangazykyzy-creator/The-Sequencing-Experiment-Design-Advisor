const fs = require('fs'); const path = require('path');
if(process.argv[2] !== '--tools') throw new Error('Expected --tools DIR');
const toolDir=path.resolve(process.argv[3]);
const createModule = require(path.join(toolDir,'minimap2.cjs'));
(async()=>{
 const out=process.argv[4]; const args=process.argv.slice(5);
 const fd=fs.openSync(out,'w'); const err=[];
 const m=await createModule({noInitialRun:true,wasmBinary:fs.readFileSync(path.join(toolDir,'minimap2.wasm')),
   getPreloadedPackage:()=>{const b=fs.readFileSync(path.join(toolDir,'minimap2.data'));return b.buffer.slice(b.byteOffset,b.byteOffset+b.byteLength)},
   locateFile:(name)=>path.join(toolDir,name),print:(line)=>fs.writeSync(fd,line+'\n'),printErr:(line)=>{err.push(line);console.error(line)}});
 m.FS.mkdir('/input');
 const mapped=[];
 for(const a of args){if(fs.existsSync(a)&&fs.statSync(a).isFile()){const n='/input/'+path.basename(a);m.FS.writeFile(n,fs.readFileSync(a));mapped.push(n);}else mapped.push(a);}
 const t=Date.now();try{m.callMain(mapped)}finally{fs.closeSync(fd);fs.writeFileSync(out+'.log',err.join('\n'));fs.writeFileSync(out+'.runtime.json',JSON.stringify({args,wasm_version:'minimap2 2.22 biowasm',elapsed_seconds:(Date.now()-t)/1000,node_rss_at_completion_bytes:process.memoryUsage().rss},null,2));}
})().catch(e=>{console.error(e);process.exit(1)});
