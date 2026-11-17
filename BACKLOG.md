# Backlog `v0.3.0` (Hallazgos no bloqueantes post-`v0.2.0`)

Lista de mejoras diferidas tras la auditoría de la ronda `7e5fc46 -> c5af733` (`v0.2.0`). Ninguna supone un riesgo de seguridad ni compromete la corrección científica actual.

## 1. Conformidad JSON-RPC / MCP
1. Devolver error JSON-RPC `-32602` (*Invalid params*) ante parámetros inválidos en lugar de respuestas de herramienta genéricas.
2. Validar estrictamente `"jsonrpc": "2.0"` y rechazar peticiones sin el campo `jsonrpc` o con `"1.0"`.
3. Ajustar el manejo de `id: null` según la especificación JSON-RPC 2.0 y no responder a notificaciones sin `id`.
4. Exigir `format` explícitamente cuando está marcado como `required` (sin *fallback* silencioso a `csv`) y declarar `outputSchema` en las definiciones de herramientas MCP.

## 2. CLI y ergonomía de sandbox
5. Validar `max_lipinski_violations` (rechazar valores negativos como `-1` en lugar de ignorarlos en silencio).
6. Evitar crear jerarquías de directorios arbitrarias al pasar `--write-root` sin aviso o confirmación explícita.
7. En modo `batch` sin `--write-root`/`--read-root`, avisar o exigir confirmación antes de escribir en rutas arbitrarias del sistema de archivos.

## 3. Informes y GUI
8. Aplicar el truncado de etiquetas largas también en el SVG BOILED-Egg (`crates/admet-report/src/svg.rs`), no solo en el lienzo de la GUI.
9. Resincronizar el sketcher 2D también tras `FilterBatch` y `ClearSession` (actualmente cubierto en `Undo`/`Redo` y selección manual).
10. Llamar a `ctx.request_repaint()` al completar exportaciones asíncronas y evitar bloquear el hilo de interfaz en importaciones pesadas.
11. Soportar historial *Undo/Redo* consistente en el editor multilínea SMILES cuando se alterna con comandos de sesión.
12. Unificar el renderizador vectorial 2D (actualmente implementado de forma independiente para GUI, SVG y PDF).
13. Comprimir streams PDF con `FlateDecode` para reducir el tamaño de lotes grandes (p. ej. ~75 MB en 5.000 moléculas).
14. Filtrar/sustituir *noncharacters* XML (`U+FFFE`, `U+FFFF`) en `escape_xml` además de los caracteres de control ASCII.
