# OS Grid Reference to WGS84 Conversion

GPT-5 vibecoded conversion from UK Ordnance Survey Grid References to and from WGS84 coordinate system.

No external script dependencies or API lookups.

Share and enjoy.

To install dependencies:

```bash
bun install
```

To run tests:

```bash
bun run test
```

To transpile TypeScript:

```bash
bun run tsc
```

To convert from OS Grid Reference to WGS84 latitude,longitude:

```
const { lat, lon } = gridRefToWgs84("TL 44982 57869");
```

To convert from WGS84 latitude,longitude to OS Grid Reference:

```
const gridRef = wgs84ToGridRef(50.805547, -0.037451998);
```
